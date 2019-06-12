#include "ORCA7.h"

#include "RooGaussian.h"
#include "RooConstVar.h"

using namespace RooFit;
using namespace O7;

//*********************************************************************************************************

// functions to  be able to use the shower energy and track direction for high-purity tracks
//---------------------------------------------------------------------------------------
Double_t O7::CustomEnergy(SummaryEvent* evt) {

  if ( evt->Get_shower_ql0() > 0.5 ) { return evt->Get_shower_energy(); }
  else                               { return evt->Get_track_energy();  }

}

TVector3 O7::CustomDir(SummaryEvent *evt) { return evt->Get_track_dir();      }
TVector3 O7::CustomPos(SummaryEvent *evt) { return evt->Get_track_pos();      }
Double_t O7::CustomBY (SummaryEvent *evt) { return evt->Get_track_bjorkeny(); }

//*********************************************************************************************************

/** Constructor */
ORCA7::ORCA7(UInt_t ResponseType) {

  // set the PID bins
  fPidBins.push_back( PidBinConf(0.0, 0.3, 0.05, 0.0, "shw") );
  fPidBins.push_back( PidBinConf(0.3, 0.7, 0.01, 0.0, "mid") );
  fPidBins.push_back( PidBinConf(0.7, 1.0+1e-5, 0.05, 0.1, "trk") );

  CreateResponses( fPidBins, ResponseType );

  // create fitutil and pdfs
  fFitUtil = new FitUtilWsyst(f_F_runtime, fResps["trk"]->GetHist3DTrue(), fResps["trk"]->GetHist3DReco(), f_F_emin, f_F_emax, 
			      f_F_ctmin, f_F_ctmax, f_F_bymin, f_F_bymax, fEffmF);
  
  for (auto P: fResps) {
    TString pdfname = "pdf_" + P.first;
    FitPDF* pdf = new FitPDF(pdfname, pdfname, fFitUtil, P.second);
    fPdfs.insert( std::make_pair(P.first, pdf) );
  }

  CreatePriors( fFitUtil );
  PrepareParameters( fFitUtil );

} // end of constructor

//*********************************************************************************************************

/** Inline function to construct the responses in the constructor*/
void ORCA7::CreateResponses(vector< O7::PidBinConf > pid_bins, UInt_t ResponseType) {

  //=========================================================
  // create a response for each PID bin
  //=========================================================
  std::map<UInt_t, TString> rmap = { {DETR, "DetResponse"}, {EVTR, "EvtResponse"}, {EVTR_EXTW1Y, "EvtResponse_externalW1Y"} };
  
  cout << "NOTICE ORCA7::CreateResponses() initialising " << rmap[ResponseType] << endl;
  
  for (auto PB: pid_bins) {

    Bool_t UseExtW1Y = ( ResponseType == EVTR_EXTW1Y );
    AbsResponse *R = NULL;

    // init shower and middle responses
    if ( PB.pid_min < 0.7 ) {

      if ( ResponseType == DETR ) {
	R = new DetResponse(DetResponse::shower, PB.name+"R",
			    f_T_ebins, f_emin, f_emax, f_T_ctbins, f_ctmin, f_ctmax, f_T_bybins, f_bymin, f_bymax,
			    f_R_ebins, f_emin, f_emax, f_R_ctbins, f_ctmin, f_ctmax, f_R_bybins, f_bymin, f_bymax);

      }
      else {
	R = new EvtResponse(EvtResponse::shower, PB.name+"R",
			    f_T_ebins, f_emin, f_emax, f_T_ctbins, f_ctmin, f_ctmax, f_T_bybins, f_bymin, f_bymax,
			    f_R_ebins, f_emin, f_emax, f_R_ctbins, f_ctmin, f_ctmax, f_R_bybins, f_bymin, f_bymax, UseExtW1Y);
      }

      R->AddCut( &SummaryEvent::Get_shower_ql0, std::greater<double>(), 0.5, true );
      
    }

    // init the track response
    else {

      if ( ResponseType == DETR ) {
	R = new DetResponse(DetResponse::customreco, PB.name+"R",
			    f_T_ebins, f_emin, f_emax, f_T_ctbins, f_ctmin, f_ctmax, f_T_bybins, f_bymin, f_bymax,
			    f_R_ebins, f_emin, f_emax, f_R_ctbins, f_ctmin, f_ctmax, f_R_bybins, f_bymin, f_bymax);
      }
      else {
	R = new EvtResponse(EvtResponse::customreco, PB.name+"R",
			    f_T_ebins, f_emin, f_emax, f_T_ctbins, f_ctmin, f_ctmax, f_T_bybins, f_bymin, f_bymax,
			    f_R_ebins, f_emin, f_emax, f_R_ctbins, f_ctmin, f_ctmax, f_R_bybins, f_bymin, f_bymax, UseExtW1Y);
      }

      R->SetObsFuncPtrs( &CustomEnergy, &CustomDir, &CustomPos, &CustomBY );
      R->AddCut( &SummaryEvent::Get_track_ql0, std::greater<double>(), 0.5, true );
      
    }

    // add common cuts
    R->AddCut(&SummaryEvent::Get_RDF_noise_score, std::less_equal<double>()   , PB.noise_cut, true);
    R->AddCut(&SummaryEvent::Get_RDF_muon_score , std::less_equal<double>()   , PB.muon_cut , true);
    R->AddCut(&SummaryEvent::Get_RDF_track_score, std::greater_equal<double>(), PB.pid_min  , true);
    R->AddCut(&SummaryEvent::Get_RDF_track_score, std::less<double>()         , PB.pid_max  , true);
    
    fResps.insert( std::make_pair(PB.name, R) );

  }

  //=========================================================
  // Fill the responses
  //=========================================================
  
  SummaryParser sp(fDataF);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i % 500000 == 0) cout << "NOTICE ORCA7::ORCA7() event " << i << " of " 
			      << sp.GetTree()->GetEntries() << endl;
    for (auto R: fResps) R.second->Fill( sp.GetEvt(i) );
  }

}

//*********************************************************************************************************

/** Inline function to construct the priors in the constructor */
void ORCA7::CreatePriors(FitUtil *F) {

  //=========================================================
  // create priors
  //=========================================================

  // need to add some priors; skew parameter priors ball-parked from Barr
  // the uncertainty there is <=20%, a prior of 20% hopefully helps to avoid the fitter going
  // to regions we know are not physical
  RooGaussian *mu_amu_prior = new RooGaussian("mu_amu_prior", "mu_amu_prior", 
					      *F->GetVar("skew_mu_amu"), RooConst(0.), RooConst(0.1) );
  RooGaussian *e_ae_prior = new RooGaussian("e_ae_prior"  , "e_ae_prior"  , 
					    *F->GetVar("skew_e_ae")  , RooConst(0.), RooConst(0.1) );
  RooGaussian *mu_e_prior = new RooGaussian("mu_e_prior"  , "mu_e_prior"  , 
					    *F->GetVar("skew_mu_e")  , RooConst(0.), RooConst(0.1) );

  // nc normalisation taken from Neutrino2018 poster
  RooGaussian *ncnorm_prior = new RooGaussian("ncnorm_prior", "ncnorm_prior", 
					      *F->GetVar("NC_norm"), RooConst(1.), RooConst(0.1) );

  fPriors.add( RooArgSet(*mu_amu_prior, *e_ae_prior, *mu_e_prior, *ncnorm_prior) );

}

//*********************************************************************************************************

/** Function to add a prior for dm31 */
void ORCA7::AddDm31Prior(Bool_t InvertedOrdering) {

  Double_t mean  = 2.525*1e-3;
  Double_t sigma = (2.622*1e-3 - 2.431*1e-3)/3;

  if ( InvertedOrdering ) {
    mean  = -2.512*1e-3 + 7.39*1e-5;
    sigma = ( (-2.413*1e-3 + 7.39*1e-5) - (-2.606*1e-3 + 7.39*1e-5) ) / 3;
  }

  // nc normalisation taken from Neutrino2018 poster
  RooGaussian *dm31_prior = new RooGaussian("dm31_prior", "dm31_prior", 
					    *fFitUtil->GetVar("Dm31"), RooConst(mean), RooConst(sigma) );
  
  fPriors.add( *dm31_prior );

}

//*********************************************************************************************************

/** Inline function to populate member vectors that differentiate between oscillation and systematic parameters */
void ORCA7::PrepareParameters(FitUtil *F) {

  // vector with pointers to oscillation parameters

  fOscPars = { F->GetVar("SinsqTh12"), F->GetVar("SinsqTh13"), F->GetVar("SinsqTh23"),
	       F->GetVar("dcp"), F->GetVar("Dm21"), F->GetVar("Dm31") };

  // populate the vector with pointers to systematic parameters

  RooArgSet parset = F->GetSet();
  TIterator *it = parset.createIterator();
  RooRealVar* var;

  while ( ( var = (RooRealVar*)it->Next() ) ) {

    if ( F->GetObs().find(var->GetName()) != NULL ) continue; // ignore observables

    if ( std::find( fOscPars.begin(), fOscPars.end(), var ) == fOscPars.end() ) {
      fSystPars.push_back( var );
    }

  }
  
  // save the systematics default values
  for (auto sv: fSystPars) {
    fSystDefault.insert( std::make_pair( sv, sv->getVal() ) );
  }

}

//*********************************************************************************************************

/** Destructor */
ORCA7::~ORCA7() {

  delete fFitUtil;
  for (auto R: fResps) if (R.second) delete R.second;
  for (auto P: fPdfs) if (P.second) delete P.second;

}

//*********************************************************************************************************

/** Set oscillation parameters to normal ordering and systematic parameters to default values*/
void ORCA7::Set_NuFit_4p0_NO() {

  FitUtil* F = fFitUtil;

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

  for (auto kv: fSystDefault) {
    F->GetVar( kv.first->GetName() )->setVal( kv.second );
  }

};

//*********************************************************************************************************

/** Set oscillation parameters to normal ordering and systematic parameters to default values*/
void ORCA7::Set_NuFit_3p2_NO() {

  FitUtil* F = fFitUtil;

  F->FreeParLims();

  F->GetVar("SinsqTh12")->setVal( 0.307 );
  F->GetVar("SinsqTh12")->setMin( 0.272 );
  F->GetVar("SinsqTh12")->setMax( 0.346 );

  F->GetVar("SinsqTh13")->setVal( 0.02206 );
  F->GetVar("SinsqTh13")->setMin( 0.01981 );
  F->GetVar("SinsqTh13")->setMax( 0.02436 );

  F->GetVar("SinsqTh23")->setVal( 0.538 );
  F->GetVar("SinsqTh23")->setMin( 0.418 );
  F->GetVar("SinsqTh23")->setMax( 0.613 );

  F->GetVar("dcp")      ->setVal( 234. * TMath::DegToRad() / TMath::Pi() );
  F->GetVar("dcp")      ->setMin( 144. * TMath::DegToRad() / TMath::Pi() );
  F->GetVar("dcp")      ->setMax( 374. * TMath::DegToRad() / TMath::Pi() );

  F->GetVar("Dm21")     ->setVal( 7.40*1e-5  );
  F->GetVar("Dm21")     ->setMin( 6.80*1e-5  );
  F->GetVar("Dm21")     ->setMax( 8.02*1e-5  );

  F->GetVar("Dm31")     ->setVal( 2.494*1e-3 );
  F->GetVar("Dm31")     ->setMin( 2.399*1e-3 );
  F->GetVar("Dm31")     ->setMax( 2.593*1e-3 );

  for (auto kv: fSystDefault) {
    F->GetVar( kv.first->GetName() )->setVal( kv.second );
  }

};

//*********************************************************************************************************

/** Set oscillation parameters to inverted ordering and systematic parameters to default values*/
void ORCA7::Set_NuFit_4p0_IO() {

  FitUtil* F = fFitUtil;

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

  for (auto kv: fSystDefault) {
    F->GetVar( kv.first->GetName() )->setVal( kv.second );
  }

}

//*********************************************************************************************************

/** Set oscillation parameters to inverted ordering and systematic parameters to default values*/
void ORCA7::Set_NuFit_3p2_IO() {

  FitUtil* F = fFitUtil;

  F->FreeParLims();

  F->GetVar("SinsqTh12")->setVal( 0.307   );
  F->GetVar("SinsqTh12")->setMin( 0.272   );
  F->GetVar("SinsqTh12")->setMax( 0.346   );

  F->GetVar("SinsqTh13")->setVal( 0.02227 );
  F->GetVar("SinsqTh13")->setMin( 0.02006 );
  F->GetVar("SinsqTh13")->setMax( 0.02452 );

  F->GetVar("SinsqTh23")->setVal( 0.554   );
  F->GetVar("SinsqTh23")->setMin( 0.435   );
  F->GetVar("SinsqTh23")->setMax( 0.616   );

  F->GetVar("dcp")      ->setVal( 278. * TMath::DegToRad() / TMath::Pi() );
  F->GetVar("dcp")      ->setMin( 192. * TMath::DegToRad() / TMath::Pi() );
  F->GetVar("dcp")      ->setMax( 354. * TMath::DegToRad() / TMath::Pi() );

  F->GetVar("Dm21")     ->setVal( 7.40*1e-5 );
  F->GetVar("Dm21")     ->setMin( 6.80*1e-5  );
  F->GetVar("Dm21")     ->setMax( 8.02*1e-5  );

  F->GetVar("Dm31")     ->setVal( -2.465*1e-3 + 7.40*1e-5 );
  F->GetVar("Dm31")     ->setMin( -2.562*1e-3 + 7.40*1e-5 );
  F->GetVar("Dm31")     ->setMax( -2.369*1e-3 + 7.40*1e-5 );

  for (auto kv: fSystDefault) {
    F->GetVar( kv.first->GetName() )->setVal( kv.second );
  }

}

//*********************************************************************************************************

void ORCA7::RandomisePars(Bool_t InvertedOrdering, Bool_t RandomiseSyst, Int_t seed) {

  TRandom rand(seed);

  //=======================================================
  // randomisation of oscillation parameters
  //=======================================================

  // set to NuFit 4.0 values and limits
  if (InvertedOrdering) Set_NuFit_4p0_IO();
  else Set_NuFit_4p0_NO();

  for (auto var: fOscPars) {
    cout << "NOTICE ORCA7::RandomisePars() randomising parameter " << var->GetName() << endl;
    Double_t mean  = var->getVal();
    Double_t sigma = ( var->getMax() - var->getMin() )/3 ;
    var->setVal( rand.Gaus( mean, sigma ) );
  }

  //=======================================================
  // randomisation of systematics
  //=======================================================

  if (!RandomiseSyst) return;

  for (auto var: fSystPars) {
    cout << "NOTICE ORCA7::RandomisePars() randomising parameter " << var->GetName() << endl;
    Double_t mean  = var->getVal();
    Double_t sigma = 0.15;
    var->setVal( rand.Gaus(mean, sigma) );
  }

}

//*********************************************************************************************************

void ORCA7::ReadFitData(TString infile) {

  TFile fin(infile, "READ");
  
  if ( !fin.IsOpen() ) {
    throw std::invalid_argument("ERROR! ORCA7::ReadFitData() cannot open file " + (string)infile);
  }

  TTree *tin = (TTree*)fin.Get("fitresult");

  if ( tin == NULL ) {
    throw std::invalid_argument("ERROR! ORCA7::ReadFitData(0 cannot find TTree fitpacket in input file");
  }

  fitpacket *packet = new fitpacket();
  tin->SetBranchAddress("fitpacket", &packet);

  for (Int_t i = 0; i < tin->GetEntries(); i++) {
    tin->GetEntry(i);
    fFPs.push_back( new fitpacket( *packet ) );
  }

  cout << "NOTICE ORCA7::ReadFitData() read " << fFPs.size() << " fitpacket's to member fFPs" << endl;

}
