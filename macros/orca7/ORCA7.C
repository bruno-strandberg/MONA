#include "ORCA7.h"

#include "RooGaussian.h"
#include "RooConstVar.h"

using namespace RooFit;
using namespace O7;

// functions to  be able to use the shower energy and track direction for high-purity tracks
//---------------------------------------------------------------------------------------
Double_t O7::CustomEnergy(SummaryEvent* evt) {

  if ( evt->Get_shower_ql0() > 0.5 ) { return evt->Get_shower_energy(); }
  else                               { return evt->Get_track_energy();  }

}

TVector3 O7::CustomDir(SummaryEvent *evt) { return evt->Get_track_dir();      }
TVector3 O7::CustomPos(SummaryEvent *evt) { return evt->Get_track_pos();      }
Double_t O7::CustomBY (SummaryEvent *evt) { return evt->Get_track_bjorkeny(); }

/** Constructor */
ORCA7::ORCA7(Bool_t ReadResponses) {

  // set the PID bins
  //-------------------------------------------------------------------------------
  fPidBins.push_back( PidBinConf(0.0, 0.3, 0.05, 0.0, "shw") );
  fPidBins.push_back( PidBinConf(0.3, 0.7, 0.01, 0.0, "mid") );
  fPidBins.push_back( PidBinConf(0.7, 1.0+1e-5, 0.05, 0.1, "trk") );

  // create a response for each PID bin
  //-------------------------------------------------------------------------------
  for (auto PB: fPidBins) {

    DetResponse *R = NULL;
      
    if ( PB.pid_min < 0.7 ) {
      R = new DetResponse(DetResponse::shower, PB.name+"R", f_R_ebins, f_R_emin, f_R_emax, 
			  f_R_ctbins, f_R_ctmin, f_R_ctmax, f_R_bybins, f_R_bymin, f_R_bymax);
      R->AddCut( &SummaryEvent::Get_shower_ql0, std::greater<double>(), 0.5, true );
    }
    else {
      R = new DetResponse(DetResponse::customreco, PB.name+"R", f_R_ebins, f_R_emin, f_R_emax, 
			  f_R_ctbins, f_R_ctmin, f_R_ctmax, f_R_bybins, f_R_bymin, f_R_bymax);
      R->SetObsFuncPtrs( &CustomEnergy, &CustomDir, &CustomPos, &CustomBY );
      R->AddCut( &SummaryEvent::Get_track_ql0, std::greater<double>(), 0.5, true );
    }

    R->AddCut(&SummaryEvent::Get_RDF_noise_score, std::less_equal<double>()   , PB.noise_cut, true);
    R->AddCut(&SummaryEvent::Get_RDF_muon_score , std::less_equal<double>()   , PB.muon_cut , true);
    R->AddCut(&SummaryEvent::Get_RDF_track_score, std::greater_equal<double>(), PB.pid_min  , true);
    R->AddCut(&SummaryEvent::Get_RDF_track_score, std::less<double>()         , PB.pid_max  , true);

    fResps.insert( std::make_pair(PB.name, R) );

  }

  // create an output name for each response
  //-------------------------------------------------------------------------------

  TString out_dir = NMHUtils::Getcwd() + "/rootfiles/"; // directory where responses are stored
  system("mkdir -p " + out_dir);

  vector< std::pair<TString, DetResponse*> > resp_names;// vector that stores the outname and the response pairs
  Bool_t ReadFromFile = ReadResponses;                  // flag to indicate whether responses should be read

  for (auto P: fResps) {
    auto R = P.second;
    TString outname = out_dir + "resp_" + R->GetRespName() + ".root";
    resp_names.push_back( std::make_pair(outname, R) );
    ReadFromFile = ReadFromFile && NMHUtils::FileExists(outname);
  }

  // read the responses from file or fill from summary data
  //-------------------------------------------------------------------------------
  if ( ReadFromFile ) {
    for (auto p: resp_names) p.second->ReadFromFile( p.first );
  }
  else {
      
    SummaryParser sp(fDataF);
    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
      if (i % 100000 == 0) cout << "NOTICE ORCA7::ORCA7() event " << i << " of " 
				<< sp.GetTree()->GetEntries() << endl;
      for (auto R: fResps) R.second->Fill( sp.GetEvt(i) );
    }
      
    for (auto p: resp_names) p.second->WriteToFile( p.first );

  }

  // create fitutil and pdfs
  //-------------------------------------------------------------------------------
  fFitUtil = new FitUtilWsyst(f_F_runtime, fResps["trk"]->GetHist3D(), f_F_emin, f_F_emax, 
			      f_F_ctmin, f_F_ctmax, f_F_bymin, f_F_bymax, fEffmF);
  
  for (auto P: fResps) {
    TString pdfname = "pdf_" + P.first;
    FitPDF* pdf = new FitPDF(pdfname, pdfname, fFitUtil, P.second);
    fPdfs.insert( std::make_pair(P.first, pdf) );
  }

  // create priors
  //-------------------------------------------------------------------------------

  // need to add some priors; skew parameter priors ball-parked from Barr
  // the uncertainty there is <=20%, a prior of 20% hopefully helps to avoid the fitter going
  // to regions we know are not physical
  RooGaussian *mu_amu_prior = new RooGaussian("mu_amu_prior", "mu_amu_prior", 
					      *fFitUtil->GetVar("skew_mu_amu"), RooConst(0.), RooConst(0.2) );
  RooGaussian *e_ae_prior = new RooGaussian("e_ae_prior"  , "e_ae_prior"  , 
					    *fFitUtil->GetVar("skew_e_ae")  , RooConst(0.), RooConst(0.2) );
  RooGaussian *mu_e_prior = new RooGaussian("mu_e_prior"  , "mu_e_prior"  , 
					    *fFitUtil->GetVar("skew_mu_e")  , RooConst(0.), RooConst(0.2) );

  // energy scale prior is a complete guess; setting 0.15, which means 2sigma is 30%, seems kind-of reasonable
  RooGaussian *escale_prior = new RooGaussian("escale_prior", "escale_prior", 
					      *fFitUtil->GetVar("E_scale"), RooConst(0.), RooConst(0.15));

  // nc normalisation taken from Neutrino2018 poster
  RooGaussian *ncnorm_prior = new RooGaussian("ncnorm_prior", "ncnorm_prior", 
					      *fFitUtil->GetVar("NC_norm"), RooConst(0.), RooConst(0.1) );

  fPriors.add( RooArgSet(*mu_amu_prior, *e_ae_prior, *mu_e_prior, *escale_prior, *ncnorm_prior) );
    
} // end of constructor

//*********************************************************************************************************

ORCA7::~ORCA7() {

  delete fFitUtil;
  for (auto R: fResps) if (R.second) delete R.second;
  for (auto P: fPdfs) if (P.second) delete P.second;

}

//*********************************************************************************************************

void ORCA7::Set_NuFit_4p0_NO(FitUtil* F) {

  F->FreeParLims();

  F->GetVar("SinsqTh12")->setVal( 0.310   );
  F->GetVar("SinsqTh12")->setMin( 0.275   );
  F->GetVar("SinsqTh12")->setMax( 0.350   );

  F->GetVar("SinsqTh13")->setVal( 0.02240 );
  F->GetVar("SinsqTh13")->setMin( 0.02044 );
  F->GetVar("SinsqTh13")->setMax( 0.02437 );

  F->GetVar("SinsqTh23")->setVal( 0.582   );
  F->GetVar("SinsqTh23")->setMin( 0.428   );
  F->GetVar("SinsqTh23")->setMax( 0.624   );

  F->GetVar("dcp")      ->setVal( 217. * TMath::DegToRad() / TMath::Pi() );
  F->GetVar("dcp")      ->setMin( 135. * TMath::DegToRad() / TMath::Pi() );
  F->GetVar("dcp")      ->setMax( 366. * TMath::DegToRad() / TMath::Pi() );

  F->GetVar("Dm21")     ->setVal( 7.39*1e-5  );
  F->GetVar("Dm21")     ->setMin( 6.79*1e-5  );
  F->GetVar("Dm21")     ->setMax( 8.01*1e-5  );

  F->GetVar("Dm31")     ->setVal( 2.525*1e-3 );
  F->GetVar("Dm31")     ->setMin( 2.431*1e-3 );
  F->GetVar("Dm31")     ->setMax( 2.622*1e-3 );

};

//*********************************************************************************************************

void ORCA7::Set_NuFit_4p0_IO(FitUtil* F) {

  F->FreeParLims();

  F->GetVar("SinsqTh12")->setVal( 0.310   );
  F->GetVar("SinsqTh12")->setMin( 0.275   );
  F->GetVar("SinsqTh12")->setMax( 0.350   );

  F->GetVar("SinsqTh13")->setVal( 0.02263 );
  F->GetVar("SinsqTh13")->setMin( 0.02067 );
  F->GetVar("SinsqTh13")->setMax( 0.02461 );

  F->GetVar("SinsqTh23")->setVal( 0.582   );
  F->GetVar("SinsqTh23")->setMin( 0.433   );
  F->GetVar("SinsqTh23")->setMax( 0.623   );

  F->GetVar("dcp")      ->setVal( 280. * TMath::DegToRad() / TMath::Pi() );
  F->GetVar("dcp")      ->setMin( 196. * TMath::DegToRad() / TMath::Pi() );
  F->GetVar("dcp")      ->setMax( 351. * TMath::DegToRad() / TMath::Pi() );

  F->GetVar("Dm21")     ->setVal( 7.39*1e-5 );
  F->GetVar("Dm21")     ->setMin( 6.79*1e-5  );
  F->GetVar("Dm21")     ->setMax( 8.01*1e-5  );

  F->GetVar("Dm31")     ->setVal( -2.512*1e-3 + 7.39*1e-5 );
  F->GetVar("Dm31")     ->setMin( -2.606*1e-3 + 7.39*1e-5 );
  F->GetVar("Dm31")     ->setMax( -2.413*1e-3 + 7.39*1e-5 );

};
