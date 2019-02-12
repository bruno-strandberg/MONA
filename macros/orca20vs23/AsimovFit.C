#include "AsimovFit.h"

#include "NMHUtils.h"
#include "FileHeader.h"
#include "SummaryParser.h"

#include "RooChi2Var.h"
#include "RooMsgService.h"

using namespace RooFit;

/** constructor */
AsimovFit::AsimovFit(Detector detector) {

  ParseInput(detector);
  InitResponses(detector);
  InitPdfs();

}

//*********************************************************************************************

/** destructor */
AsimovFit::~AsimovFit() {

  if (fTrkResp) delete fTrkResp;
  if (fShwResp) delete fShwResp;

  if (fCategs)  delete fCategs;
  if (fFitUtil) delete fFitUtil;
  if (fTrkPdf)  delete fTrkPdf;
  if (fShwPdf)  delete fShwPdf;
  if (fSimPdf)  delete fSimPdf;
  if (fFitPdf)  delete fFitPdf;
  if (fTh13C)   delete fTh13C;

}

//*********************************************************************************************

/** Function to set simulation and effective mass files, depending on detector type.
    \param detector  uint that specifies the detector
*/
void AsimovFit::ParseInput(Detector detector) {

  fDetStr   = fDetStrings[detector];
  fDataFile = fSimFiles[detector];
  fEffmFile = fEffmFiles[detector];

  // check that datatag's match in the simulation file and the effective mass file
  FileHeader h1("h1");
  FileHeader h2("h2");

  h1.ReadHeader(fDataFile);
  h2.ReadHeader(fEffmFile);
  TString datatag = "datatag";

  if ( h1.GetParameter(datatag) != h2.GetParameter(datatag) ) {
    throw std::invalid_argument("ERROR! AsimovFit::ParseInput() data-tag mismatch, found tags: " + 
				(string)h1.GetParameter(datatag) + " and " + (string)h2.GetParameter(datatag) );
  }

}

//*********************************************************************************************

/** Function that initialises and fills responses for tracks and showers, depending on the detector type.
    \param detector  Enum of type `Detector` that specifies detector type
*/
void AsimovFit::InitResponses(Detector detector) {

  TString trkresp = fDetStr + "_trkresp";
  TString shwresp = fDetStr + "_shwresp";

  // track response - same irrespective of the detector
  fTrkResp = new DetResponse(DetResponse::track, trkresp, fNebins, 1, 100, fNctbins, -1, 1, fNbybins, 0, 1);
  fTrkResp->AddCut( &SummaryEvent::Get_track_ql1      , std::greater<double>(), 0.5, true );
  fTrkResp->AddCut( &SummaryEvent::Get_RDF_track_score, std::greater<double>(), 0.6, true );

  // shower response - the quality level cut is different for different geometries
  // due to bug-matching, see `apps/data_sorting/PIDGammaToSummary.C`
  fShwResp = new DetResponse(DetResponse::shower, shwresp, fNebins, 1, 100, fNctbins, -1, 1, fNbybins, 0, 1);
  fShwResp->AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), 0.6, true);
  
  switch (detector) {
   
  case ORCA20:
    fShwResp->AddCut( &SummaryEvent::Get_shower_ql2, std::greater<double>(), 0.5, true);
    break;
  case ORCA23:
    fShwResp->AddCut( &SummaryEvent::Get_shower_ql1, std::greater<double>(), 0.5, true);
    break;
    
  default:
    throw std::invalid_argument("ERROR! AsimovFit::InitResponses() unknown detector."); 
  }
  
  //-----------------------------------------------------------------------
  // fill; if possible read from file
  //-----------------------------------------------------------------------
  TString dir    = NMHUtils::Getcwd() + "/resps/";
  TString syscmd = "mkdir -p " + dir;

  if ( system(syscmd) != 0 ) { 
    cout << "NOTICE AsimovFit::InitResponses() system() returned non-zero to cmd " << syscmd << endl;
  }

  TString trk_name = dir + fTrkResp->Get_RespName() + ".root";
  TString shw_name = dir + fShwResp->Get_RespName() + ".root";  

  if ( !NMHUtils::FileExists(trk_name) || ! NMHUtils::FileExists(shw_name) ) {

    SummaryParser sp( fDataFile );

    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {

      if ( i % 1000000 == 0) cout << "Event: " << i << " of " << sp.GetTree()->GetEntries() << endl;

      fTrkResp->Fill( sp.GetEvt(i) );
      fShwResp->Fill( sp.GetEvt(i) );

    }

    fTrkResp->WriteToFile( trk_name );
    fShwResp->WriteToFile( shw_name );

  }
  else {

    fTrkResp->ReadFromFile( trk_name );
    fShwResp->ReadFromFile( shw_name );    

  }
  cout << "NOTICE AsimovFit::InitResponses() responses ready" << endl;

}

//*********************************************************************************************

/** Function to initialise PDFs */
void AsimovFit::InitPdfs() {

  using namespace RooFit;

  if ( fTrkResp == NULL || fShwResp == NULL ) {
    throw std::logic_error("ERROR! AsimovFit::InitPdfs() called before AsimovFit::InitResponses()");
  }

  // init fit utility and pdf's for tracks and showers
  fFitUtil = new FitUtil( fOpTime, fTrkResp->GetHist3D(), fFitEmin, fFitEmax, 
			  fFitCTmin, fFitCTmax, fFitBYmin, fFitBYmax, fEffmFile);

  fTrkPdf = new FitPDF(fDetStr + "_trkpdf", fDetStr + "_trkpdf", fFitUtil, fTrkResp);
  fShwPdf = new FitPDF(fDetStr + "_shwpdf", fDetStr + "_shwpdf", fFitUtil, fShwResp);

  // create a pdf for fitting tracks and showers simultaneously
  fCategs = new RooCategory("categs", "Data categories");
  fCategs->defineType(fTrkCateg);
  fCategs->defineType(fShwCateg);

  fSimPdf = new RooSimultaneous(fDetStr + "_simpdf", fDetStr + "_simpdf", *fCategs);
  fSimPdf->addPdf( *fTrkPdf, fTrkCateg );
  fSimPdf->addPdf( *fShwPdf, fShwCateg );

  // create a constraint on sinsqth13 and a pdf that includes the constrain
  fTh13C = new RooGaussian("Th13C","Th13C", *(fFitUtil->GetVar("SinsqTh13")),RooConst(fTh13mean),RooConst(fTh13sigma));

  // create a 'fitpdf', this is the simultaneous pdf that has an additional constraint on th13
  fFitPdf = new RooProdPdf(fDetStr + "_fitpdf", fDetStr + "_fitpdf", RooArgSet( *fSimPdf, *fTh13C) ) ;
  
}

//*********************************************************************************************

/** Function to create expectation values at different oscillation parameter values.
    theta-23 is the first because often the other's run at default values.
    \param th23   sin^2(th23) value
    \param th12   sin^2(th12) value
    \param th23   sin^2(th23) value
    \param dcp    delta-CP value, (0--2)
    \param dm21   delta-m21^2 value
    \param dm31   delta-m31^2 value
    \return       fitpacket with members pointing to created data
*/
fitpacket AsimovFit::CreateData(Bool_t NOdata, Double_t th23, Double_t th12, Double_t th13, Double_t dcp, Double_t dm21, Double_t dm31) {

  // set parameters
  //----------------

  if ( !NOdata ) fFitUtil->SetIOcentvals();

  if (th23 != fNOTSET) fFitUtil->GetVar("SinsqTh23")->setVal(th23);
  if (th12 != fNOTSET) fFitUtil->GetVar("SinsqTh12")->setVal(th12);
  if (th13 != fNOTSET) fFitUtil->GetVar("SinsqTh13")->setVal(th13);
  if (dcp  != fNOTSET) fFitUtil->GetVar("dcp")->setVal(dcp);
  if (dm21 != fNOTSET) fFitUtil->GetVar("Dm21")->setVal(dm21);
  if (dm31 != fNOTSET) { 
    
    if ( (NOdata && dm31 < 0 ) || (!NOdata && dm31 > 0 ) ) {
      throw std::invalid_argument("ERROR! AsimovFit::CreateData() The NOdata flag and dm31 input value conflict.");
    }
    
    fFitUtil->GetVar("Dm31")->setVal(dm31);
  }

  // create histograms and a combined dataset
  //-----------------------------------------

  TString trkname = fDetStr + "_trk";
  TString shwname = fDetStr + "_shw";
  TH3D* trkH = fTrkPdf->GetExpValHist();
  TH3D* shwH = fShwPdf->GetExpValHist();
  trkH->SetNameTitle(trkname, trkname);
  shwH->SetNameTitle(shwname, shwname);
  trkH->SetDirectory(0);
  shwH->SetDirectory(0);

  RooDataHist *combd = CombineData( trkH, shwH );
    
  // snapshot of the oscillation parameters
  //----------------------------------------
  RooArgSet *pars = (RooArgSet*)fFitUtil->GetSet().snapshot(kTRUE);

  // set pointers for the fitpacket
  //-------------------------------
  fitpacket fp;
  fp.fDetString = fDetStr;
  fp.fTrkH      = trkH;
  fp.fShwH      = shwH;
  fp.fCombH     = combd;
  fp.fParData   = pars;

  return fp;
}

//*********************************************************************************************

/** Function to combine the track and shower histograms to a combined dataset in RooFit 
    \param trkH Pointer to track histogram
    \param shwH Pointer to shower histogram
    \return A RooDataHist with combined data
 */
RooDataHist* AsimovFit::CombineData(TH1* trkH, TH1* shwH) {

  // category <--> data map
  std::map<string, TH1* > hist_map = { { (string)fTrkCateg, trkH }, 
				       { (string)fShwCateg, shwH } };    
  // create combined data sets
  TString combname = fDetStr + "_comb";
  RooMsgService::instance().setGlobalKillBelow(FATAL); // don't want to see messages at read-in
  RooDataHist *combd = new RooDataHist(combname, combname, fFitUtil->GetObs(), *fCategs, hist_map);
  RooMsgService::instance().setGlobalKillBelow(INFO);

  return combd;

}

//*********************************************************************************************

/** Function that fits the created data with inverted ordering in theta-23 two quadrants.
    \param fp   fitpacket returned by AsimovFit::CreateData()
 */
void AsimovFit::FitData(fitpacket& fp) {

  // set up a pointer to a function to call resetting the parameters to central values

  std::function< void(FitUtil&) > ResetToCentral;
  RooRealVar *dm31 = (RooRealVar*)fp.fParData->find("Dm31");
  if ( !dm31 ) throw std::logic_error("ERROR! AsimovFit::FitData() cannot find dm31");

  if ( dm31->getVal() > 0 ) { ResetToCentral = &FitUtil::SetIOcentvals; }
  else                      { ResetToCentral = &FitUtil::SetNOcentvals; }

  // fix solar parameters
  fFitUtil->GetVar("Dm21")->setConstant(kTRUE);
  fFitUtil->GetVar("SinsqTh12")->setConstant(kTRUE);

  // create quantiles for theta23
  fFitUtil->GetVar("SinsqTh23")->setRange("firstq" , 0., 0.5);
  fFitUtil->GetVar("SinsqTh23")->setRange("secondq", 0.5, 1.);

  RooDataHist *d = fp.fCombH;
  if (d == NULL) {
    throw std::invalid_argument("ERROR! AsimovFit::FitData() no combined data in the input fit-packet. The input fitpacket should be created with AsimovFit::CreateData()");
  }

  ResetToCentral(*fFitUtil);
  fFitUtil->GetVar("SinsqTh23")->setVal(0.4);
  RooFitResult *res_1q = fFitPdf->chi2FitTo( *d, Save(), Range("firstq"), DataError(RooAbsData::Expected) );

  ResetToCentral(*fFitUtil);
  fFitUtil->GetVar("SinsqTh23")->setVal(0.6);
  RooFitResult *res_2q = fFitPdf->chi2FitTo( *d, Save(), Range("secondq"), DataError(RooAbsData::Expected) );

  fp.fRes_1q = res_1q;
  fp.fRes_2q = res_2q;

}

//*********************************************************************************************

/** Function to write a vector of fit packets to ROOT TTree 
    \param outfile  Name of the output file
    \param fps      Vector of fit packets
 */
void AsimovFit::WriteToFile(TString outfile, vector<fitpacket> fps) {

  TFile fout(outfile, "RECREATE");
  TTree tout("fitresult", "TTree of fitpacket structures");
  
  fitpacket *packet = new fitpacket();
  tout.Branch( "fitpacket", &packet );

  for (auto fp: fps) {
    packet = &fp;
    tout.Fill();
  }
  
  tout.Write();
  fout.Close();
  
}

//*********************************************************************************************

/** Function to read fitpackets from a ROOT file.
    
    Only fit packets that match the detector string of this instance are read to memory for further manipulation.

    \param infile File where the fitpackets are stored
 */
void AsimovFit::ReadFromFile(TString infile) {

  TFile *fin = new TFile(infile, "READ");
  
  if ( !fin->IsOpen() ) {
    throw std::invalid_argument("ERROR! AsimovFit::ReadFromFile() could not open file " + (string)infile);
  }

  TTree *tin = (TTree*)fin->Get("fitresult");

  if (tin == NULL) {
    throw std::invalid_argument("ERROR! AsimovFit::ReadFromFile() could not find TTree fitresult");
  }

  fitpacket *packet = new fitpacket();
  tin->SetBranchAddress("fitpacket", &packet);

  Int_t ignored = 0;
  for (Int_t i = 0; i < tin->GetEntries(); i++) {
    tin->GetEntry(i);

    if ( packet->fDetString == fDetStr ) {

      fFPs.push_back( new fitpacket( *packet ) );

      // the combined dataset needs to be re-created because otherwise it is not related to the
      // fCategs of this class instance
      fFPs.back()->fCombH = CombineData( fFPs.back()->fTrkH, fFPs.back()->fShwH );
    }
    else { ignored++; }

  }

  cout << "NOTICE AsimovFit::ReadFromFile() Read in: " << fFPs.size() << " fit packets out of "
       << tin->GetEntries() << ", " << ignored << " belonged to a different detector config." << endl;

  fin->Close();
  if (fin) delete fin;

}

//*********************************************************************************************

/** Function to scan the likelihood in a given variable in a given range
    \param var_name     Name of the variable as defined in FitUtil
    \param var_min      Start value of the scan
    \param var_max      Stop value of the scan
    \param th23         sin^(th23) value for the data with respect to which the scan is performed
 */
void AsimovFit::LikelihoodScan(TString var_name, Double_t var_min, Double_t var_max, Double_t th23) {

  auto data = CreateData( th23 );
 
  // set to IO central values, fix solar parameters
  fFitUtil->SetIOcentvals();
  fFitUtil->GetVar("Dm21")->setConstant(kTRUE);
  fFitUtil->GetVar("SinsqTh12")->setConstant(kTRUE);

  // create LLH calculator
  RooNLLVar nll( "nll", "nll", *fFitPdf, *(data.fCombH), NumCPU(4) );
  TString title = "LLH scan wrt " + (TString)fFitUtil->GetVar(var_name)->GetName();
  RooPlot *frame = fFitUtil->GetVar(var_name)->frame( Range(var_min, var_max), Title(title) );
  nll.plotOn( frame, LineColor(kRed), ShiftToZero() );
  
  new TCanvas("c1","c1",1);
  frame->Draw();

}

//*********************************************************************************************

/** Function to draw a contour plot at a certain th23 value in two variables
    \param th23_data    th23 value at which the data is generated
    \param var1         FitUtil parameter name 1, default SinsqTh23
    \param var2         FitUtil parameter name 2, default SinsqTh13
*/
void AsimovFit::Contour(Double_t th23_data, TString var1, TString var2) {
  
  auto data = CreateData(th23_data, th23_data, 0.1);

  // set to IO central values, fix solar parameters
  fFitUtil->SetIOcentvals();
  fFitUtil->GetVar("Dm21")->setConstant(kTRUE);
  fFitUtil->GetVar("SinsqTh12")->setConstant(kTRUE);

  // create LLH calculator
  RooNLLVar nll( "nll", "nll", *fFitPdf, *(data.fCombH), NumCPU(4) );
  RooMinuit min( nll );
  RooPlot *contour = min.contour( *(fFitUtil->GetVar(var1)), *(fFitUtil->GetVar(var2)), 2 );

  new TCanvas("ccont","ccont",1);
  contour->Draw();
    
}

//*********************************************************************************************

/** Function to calculate the chi2 values between the model at the fit result and the data as stored in the fitpacket.

    Note on errors and their effect on chi2: RooFit treats the data with errors "RooAbsData::Expected", which, as much as I undertand, means that it assumes the error in the data bin to be sqrt(bin_content). The bin content is the expectation value in a given ordering, the selected error aims to mimic what we do in NMHUtils::Asymmetry, and does not use the MC statistical error that is stored in SumW2 of the expectation value histograms created with FitPDF::GetExpValHist(). The model does not have an associated uncertainty. I tested some of the steps against TH1::Chi2Test function, which seemed to suggest this is a sound procedure.

    \param fp Fit packet with data, parameters at data and parameters at fit
    \param q1 Boolean to select either first (q1=true) or second (q2=false) th23 quadrant
    \return a tuple with: <0> combined chi2, including constraint on th13 ; <1> chi2 on tracks ; <2> chi2 on showers. <1-2> do not include info on th13 constraint.
*/
std::tuple<Double_t, Double_t, Double_t> AsimovFit::GetChi2(fitpacket &fp, Bool_t q1) {

  // check that fitpacket was created with the same detector selection
  if ( fDetStr != fp.fDetString ) {
    throw std::invalid_argument("ERROR! AsimovFit::GetChi2() detector string mismatch");
  }

  // prepare data
  RooMsgService::instance().setGlobalKillBelow(FATAL); // don't want to see messages...
  RooDataHist trkd_rf("trkd_rf","trkd_rf", fFitUtil->GetObs(), Import(*(fp.fTrkH)) );
  RooDataHist shwd_rf("shwd_rf","shwd_rf", fFitUtil->GetObs(), Import(*(fp.fShwH)) );
  RooMsgService::instance().setGlobalKillBelow(INFO); // want to see them again...

  //select the fit result depending on the quarter and set model to fit result values
  RooFitResult *fitres;
  if (q1) fitres = fp.fRes_1q;
  else fitres = fp.fRes_2q;

  SetModelToFitresult(fitres);

  // create chi2 variables
  RooChi2Var chi2_trk("chi2_trk","chi2_trk", *fTrkPdf, trkd_rf, DataError(RooAbsData::Expected) );
  RooChi2Var chi2_shw("chi2_shw","chi2_shw", *fShwPdf, shwd_rf, DataError(RooAbsData::Expected) );

  // additional chi2 term from the constraint on theta-13
  Double_t th13_fitted = ( (RooRealVar*)fitres->floatParsFinal().find("SinsqTh13") )->getVal();
  Double_t chi2_th13 = TMath::Power( (th13_fitted - fTh13mean)/fTh13sigma, 2 );

  Double_t trkX2  = chi2_trk.getVal();
  Double_t shwX2  = chi2_shw.getVal();
  Double_t combX2 = trkX2 + shwX2 + chi2_th13;

  return std::make_tuple( combX2, trkX2, shwX2 );

}

//*********************************************************************************************

std::tuple<Double_t, Double_t, Double_t> AsimovFit::GetAsym(fitpacket &fp, Bool_t q1) {

  //prep pdf's
  if (q1) SetModelToFitresult(fp.fRes_1q);
  else SetModelToFitresult(fp.fRes_2q);

  // get the histograms of other ordering as established by the minimizer
  TH3D* trk_fit = fTrkPdf->GetExpValHist();
  TH3D* shw_fit = fShwPdf->GetExpValHist();

  // calculate asymmetries
  auto A_trk = NMHUtils::Asymmetry(fp.fTrkH, trk_fit, "asym_trk");
  auto A_shw = NMHUtils::Asymmetry(fp.fShwH, shw_fit, "asym_shw");

  Double_t combA = std::get<1>(A_trk) * std::get<1>(A_trk) + std::get<1>(A_shw) * std::get<1>(A_shw);

  return std::make_tuple(combA, TMath::Power( std::get<1>(A_trk), 2 ), TMath::Power( std::get<1>(A_shw), 2 ) );

}

//*********************************************************************************************

/** Set model to results from the fit.*/
void AsimovFit::SetModelToFitresult(RooFitResult *fr) {

  TIterator* constit = fr->constPars().createIterator();
  TIterator* floatit = fr->floatParsFinal().createIterator();

  RooRealVar *constpar, *floatpar;

  while ( ( constpar = (RooRealVar*)constit->Next() ) ) {
    //cout << "NOTICE: setting constpar " << constpar->GetName() << " to " << constpar->getVal() << endl;
    fFitUtil->GetVar( constpar->GetName() )->setVal( constpar->getVal() );
  }

  while ( ( floatpar = (RooRealVar*)floatit->Next() ) ) {
    //cout << "NOTICE: setting floatpar " << floatpar->GetName() << " to " << floatpar->getVal() << endl;
    fFitUtil->GetVar( floatpar->GetName() )->setVal( floatpar->getVal() );
  }

}

//*********************************************************************************************

/** Function to find a fit packet from member fFPs vector */
fitpacket* AsimovFit::FindPacket(Double_t sinsqth23, Detector det) {

  TString detstr = fDetStrings[det];
  
  fitpacket *fp = NULL;

  for (auto f: fFPs) {

    Double_t _sinsqth23 = ( (RooRealVar*)f->fParData->find("SinsqTh23") )->getVal();

    if ( f->fDetString == detstr && _sinsqth23 == sinsqth23 ) {
      fp = f;
      break;
    }

  }

  if (fp == NULL) {
    cout << "WARNING! AsimovFit::FindPacket could not find the specified fit result and returns NULL" << endl;
  }

  return fp;
}
