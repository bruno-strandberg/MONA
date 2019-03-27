#include "ORCA7.h"

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
    TString outname = out_dir + "resp_" + R->Get_RespName() + ".root";
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


} // end of constructor

//*********************************************************************************************************

ORCA7::~ORCA7() {

  delete fFitUtil;
  for (auto R: fResps) if (R.second) delete R.second;
  for (auto P: fPdfs) if (P.second) delete P.second;

}
