#include "ORCA7.h"

#include "NMHUtils.h"
#include "FitUtilWsyst.h"
#include "FitPDF.h"

#include "RooRandom.h"
#include "RooArgSet.h"
#include "RooSimultaneous.h"
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooNLLVar.h"
#include "RooDataHist.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "RooProdPdf.h"

#include "TFile.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "Math/QuantFunc.h"

#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

using namespace RooFit;

int main(const int argc, const char** argv) {

  Bool_t   InvertedOrdering;
  Bool_t   IncludeSystematics;
  Int_t    seed;
  Int_t    ncpu;
  string   outfile;
  Double_t sinsqth23;
  Bool_t   addDm31Prior;
  
  try {

    JParser<> zap("This application creates a pseudo-experiment");

    zap['i'] = make_field(InvertedOrdering, "Inverted NMO");
    zap['s'] = make_field(IncludeSystematics, "Include systematic parameters with priors");
    zap['S'] = make_field(seed, "Seed for pseudo-experiments") = 416;
    zap['N'] = make_field(ncpu, "Number of CPUs for minimisation") = 1;
    zap['o'] = make_field(outfile, "Output file") = "fit.root";
    zap['t'] = make_field(sinsqth23, "SinsqTh23 value at which the 'data' is created") = 0.58;
    zap['p'] = make_field(addDm31Prior, "Add external prior for dm31 from NuFit");

    if ( zap.read(argc, argv) != 0 ) return 1;
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  //=====================================================================================
  // initialisation
  //=====================================================================================

  TStopwatch timer;

  // init the class with binning info, PID ranges and responses
  ORCA7 o7( kTRUE );
  FitUtilWsyst *fu = o7.fFitUtil;
  auto pdfs = o7.fPdfs;

  //=====================================================================================
  // create data, fix some of the parameters
  //=====================================================================================

  if (InvertedOrdering) { o7.Set_NuFit_4p0_IO(); }
  else                  { o7.Set_NuFit_4p0_NO(); }

  fu->GetVar("SinsqTh23")->setVal( sinsqth23 ); // theta-23 value for data creation
    
  std::map< string, TH1* > exps;     // map with response name <--> pseudoexperiment histogram
  std::map< string, TH1* > expected; // map with response name <--> expectation value histogram

  for (auto P: pdfs) {
    FitPDF *pdf = P.second;
    pdf->SetSeed( seed );

    TString hname1 = "exp_" + P.first;
    TString hname2 = "expected_" + P.first;
    TH3D   *exp    = pdf->SimplePseudoExp( hname1, kFALSE );
    TH3D   *expct  = pdf->GetExpValHist(); 
    expct->SetNameTitle(hname2, hname2);

    exps.insert( std::make_pair( (string)P.first, (TH1*)exp ) );
    expected.insert( std::make_pair( (string)P.first, (TH1*)expct ) );
  }

  RooArgSet *pars = (RooArgSet*)fu->GetSet().snapshot(kTRUE); // save the par values at which data is created
  
  if (!IncludeSystematics) { // fix the systematic parameters if they are not included
    for (auto var: o7.fSystPars) {
      var->setConstant(kTRUE);
    }
  }

  // fix tau normalisation in this analysis
  fu->GetVar("Tau_norm")->setConstant(kTRUE);

  // always fix these osc pars, no sensitivity
  fu->GetVar("SinsqTh12")->setConstant(kTRUE);
  fu->GetVar("SinsqTh13")->setConstant(kTRUE);
  fu->GetVar("dcp")->setConstant(kTRUE);
  fu->GetVar("Dm21")->setConstant(kTRUE);
  
  //=====================================================================================
  // configure simultaneous fitting
  //=====================================================================================

  RooCategory categ("categ","data categories");
  for (auto &E: exps) { categ.defineType( (TString)E.first ); }

  RooDataHist comb("comb","combined data", fu->GetObs(), categ, exps);  
  RooDataHist combExpected("combExpected","combined expectation data", fu->GetObs(), categ, expected);  

  RooSimultaneous simPdf("simPdf","simultaneous pdf", categ);
  for (auto &P: pdfs) { simPdf.addPdf( *P.second, P.first ); }

  if ( addDm31Prior ) o7.AddDm31Prior( InvertedOrdering );

  //=====================================================================================
  // fit in first and second quadrant
  //=====================================================================================

  // for fit start, set the osc parameters to NuFit central values
  if (InvertedOrdering) { o7.Set_NuFit_4p0_IO(); }
  else                  { o7.Set_NuFit_4p0_NO(); }
  fu->FreeParLims();
  fu->GetVar("SinsqTh23")->setVal(0.45);
  fu->GetVar("SinsqTh23")->setMin(0.0);
  fu->GetVar("SinsqTh23")->setMax(0.5);

  RooFitResult *fitres_q1 = simPdf.fitTo( comb, Save(kTRUE), NumCPU(ncpu), ExternalConstraints( o7.fPriors ) );

  // reset the parameters to the same starting point,except for theta-23
  if (InvertedOrdering) { o7.Set_NuFit_4p0_IO(); }
  else                  { o7.Set_NuFit_4p0_NO(); }
  fu->FreeParLims(); // release the limits
  fu->GetVar("SinsqTh23")->setVal(0.55);
  fu->GetVar("SinsqTh23")->setMin(0.5);
  fu->GetVar("SinsqTh23")->setMax(1.0);

  RooFitResult *fitres_q2 = simPdf.fitTo( comb, Save(kTRUE), NumCPU(ncpu), ExternalConstraints( o7.fPriors ) );

  //=====================================================================================
  // print results
  //=====================================================================================

  RooFitResult* result;
  (fitres_q1->minNll() < fitres_q2->minNll() ) ? result = fitres_q1 : result = fitres_q2;

  TIterator *floatit = result->floatParsFinal().createIterator();
  TIterator *constit = result->constPars().createIterator();
  RooRealVar *floatpar, *constpar;

  while ( ( floatpar = (RooRealVar*)floatit->Next() ) ) {
    RooRealVar *init = (RooRealVar*)pars->find( floatpar->GetName() );
    cout << "NOTICE fit: fitted parameter: " << floatpar->GetName() << ", true and fitted: " 
	 << init->getVal() << "\t" << floatpar->getVal() << endl;
  }

  while ( ( constpar = (RooRealVar*)constit->Next() ) ) {
    cout << "NOTICE fit: fixed parameter: " << constpar->GetName() << endl;
  }  
    
  cout << "Total run time: " << (Double_t)timer.RealTime() << endl;

  //=====================================================================================
  // scan the likelihood both with respect to the expectation data and pseudo-experiment data
  //=====================================================================================

  RooRealVar* V = fu->GetVar("SinsqTh23");
  RooPlot frame( *V, 0.3, 0.7, 40 );

  RooAbsReal* nll1 = simPdf.createNLL( comb, ExternalConstraints( o7.fPriors ), NumCPU(ncpu) );
  nll1->plotOn( &frame, ShiftToZero(), Name("NLLexp"), LineColor(kRed) );

  RooAbsReal* nll2 = simPdf.createNLL( combExpected, ExternalConstraints( o7.fPriors ), NumCPU(ncpu) );
  nll2->plotOn( &frame, ShiftToZero(), Name("NLLexpct"), LineColor(kBlue) );

  TGraph *gnll_1 = (TGraph*)frame.findObject( "NLLexp" )->Clone();
  TGraph *gnll_2 = (TGraph*)frame.findObject( "NLLexpct" )->Clone();
  TString prefix = "-log(L) scan in sin^{2}#theta_{23}";
  gnll_1->SetNameTitle( prefix + ", pseudo-exp", prefix + ", pseudo-exp" );
  gnll_2->SetNameTitle( prefix + ", expectation", prefix + ", expectation" );

  vector< TGraph* > graphs = {gnll_1, gnll_2};
  for (auto g: graphs) {
    g->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
    g->GetYaxis()->SetTitle("-log(L)");
  }

  //=====================================================================================
  // write to output a TTree of fitpackets. Although there is only one entry, this simplifies
  // the following analysis as one can hadd 1k fit outputs together and then loop over the tree
  //=====================================================================================

  TFile fout((TString)outfile, "RECREATE");
  TTree tout("fitresult", "TTree of fitpackets");

  fitpacket *fp = new fitpacket();
  tout.Branch("fitpacket", &fp);

  fp->fShw     = (TH3D*)exps["shw"];
  fp->fMid     = (TH3D*)exps["mid"];
  fp->fTrk     = (TH3D*)exps["trk"];
  fp->fParData = pars;
  fp->fRes_1q  = fitres_q1;
  fp->fRes_2q  = fitres_q2;
  fp->fLLHscan_exp      = gnll_1;
  fp->fLLHscan_expected = gnll_2;
  fp->fSeed    = seed;

  tout.Fill();
  tout.Write();
  fout.Close();

}
