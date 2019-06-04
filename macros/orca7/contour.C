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
  Double_t sinsqth23;
  Bool_t   IncludeSystematics;
  Bool_t   IncludePriors;
  Int_t    seed;
  Int_t    ncpu;
  string   outfile;
  Bool_t   evtResp;
  Bool_t   llhScans;
  Bool_t   extendedFit;
  Bool_t   nufit3p2;
  
  try {

    JParser<> zap("This application creates contour plots for ORCA7");

    zap['i'] = make_field(InvertedOrdering, "Inverted NMO");
    zap['t'] = make_field(sinsqth23, "SinsqTh23 value at which the 'data' is created. NuFit4.0 = 0.58, NuFit3.2 = 0.538") = 0.58;
    zap['s'] = make_field(IncludeSystematics, "Include systematic parameters");
    zap['p'] = make_field(IncludePriors, "Include priors for systematics");
    zap['S'] = make_field(seed, "Seed for parameter randomisation") = 416;
    zap['N'] = make_field(ncpu, "Number of CPUs for minimisation") = 1;
    zap['o'] = make_field(outfile, "Output file") = "contour.root";
    zap['e'] = make_field(evtResp, "Use event-by-event detector response");
    zap['l'] = make_field(llhScans, "In addition to the contour, perform LLH scans in SinsqTh23");
    zap['E'] = make_field(extendedFit, "Perform and extended likelihood fit, in which case also the overall normalisation is included in the LLH");
    zap['x'] = make_field(nufit3p2, "Use NuFit3.2 values, instead of NuFit4.0. Note that th23 needs to be configured separately on the command line");

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
  ORCA7 o7( kTRUE, evtResp );
  FitUtilWsyst *fu = o7.fFitUtil;
  auto pdfs = o7.fPdfs;

  //=====================================================================================
  // create data
  //=====================================================================================

  // set the osc parameters to NuFit central values and systematics to defaults
  if (nufit3p2) {

    if (InvertedOrdering) { o7.Set_NuFit_3p2_IO(); }
    else                  { o7.Set_NuFit_3p2_NO(); }

  }
  else {

    if (InvertedOrdering) { o7.Set_NuFit_4p0_IO(); }
    else                  { o7.Set_NuFit_4p0_NO(); }

  }

  // theta-23 is set to the input value
  fu->GetVar("SinsqTh23")->setVal( sinsqth23 );
    
  // map with response name <--> expectation value histogram
  std::map< string, TH1* > exps;

  for (auto P: pdfs) {
    TString hname = "exp_" + P.first;
    FitPDF *pdf = P.second;
    TH3D   *exp = pdf->GetExpValHist();
    exp->SetNameTitle( hname, hname );
    exps.insert( std::make_pair( (string)P.first, (TH1*)exp ) );

  }

  // save the parameter values at which data is created
  RooArgSet *pars = (RooArgSet*)fu->GetSet().snapshot(kTRUE);
  
  //=====================================================================================
  // manipulate parameters for fit start
  //=====================================================================================

  fu->FreeParLims(); // release the limits
  
  // fix the systematic parameters if they are not included
  if (!IncludeSystematics) {
    for (auto var: o7.fSystPars) {
      var->setConstant(kTRUE);
    }
  }

  // tau normalisation is not included in the analysis, fix
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

  RooSimultaneous simPdf("simPdf","simultaneous pdf", categ);
  for (auto &P: pdfs) { simPdf.addPdf( *P.second, P.first ); }

  //=====================================================================================
  // create likelihood scanners - individual for each PID range and simultaneous
  //=====================================================================================

  RooRealVar* V = fu->GetVar("SinsqTh23");
  RooPlot frame( *V, 0.3, 0.7, 40 );
  frame.SetNameTitle( "-log(L) scan in sin^{2}#theta_{23}", "-log(L) scan in sin^{2}#theta_{23}" );
  frame.GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  frame.GetYaxis()->SetTitle("-log(L)");
  
  vector< RooDataHist* > rfdata;         // vector to store pointers to RooFit data, for clean-up
  std::map< TString, RooAbsReal* > nlls; // map with response name <--> nll variable

  for (auto &P: pdfs) {

    TH3D   *hin = (TH3D*)exps[ (string)P.first ];
    FitPDF *pdf = P.second;

    if ( extendedFit ) pdf->IncludeNorm();
    
    TString hname = TString("rf_") + (TString)hin->GetName();
    RooDataHist *rfh = new RooDataHist( hname, hname, fu->GetObs(), Import(*hin) );
    
    RooAbsReal *nll;
    if (IncludePriors) { nll = pdf->createNLL( *rfh, ExternalConstraints( o7.fPriors ), NumCPU(ncpu) ); }
    else { nll = pdf->createNLL( *rfh, NumCPU(ncpu) ); }

    nlls.insert( std::make_pair( P.first, nll )  );
    rfdata.push_back( rfh );
  }

  // create a LLH scanner from simultaneous fitting, add to map
  RooAbsReal *nll_sim;
  if (IncludePriors) { nll_sim = simPdf.createNLL( comb, ExternalConstraints( o7.fPriors ), NumCPU(ncpu) ); }
  else { nll_sim = simPdf.createNLL( comb, NumCPU(ncpu) ); }
  
  nlls.insert( std::make_pair("sim", nll_sim) );

  //=====================================================================================
  // Minimize and create contour plot
  //=====================================================================================

  // need to offset these for the fitter to map out the surroundings of the minimum
  fu->GetVar("SinsqTh23")->setVal( fu->GetVar("SinsqTh23")->getVal() + gRandom->Uniform(-0.02, 0.02) );
  fu->GetVar("Dm31")->setVal( fu->GetVar("Dm31")->getVal() + gRandom->Uniform(-0.00003, 0.00003) );

  RooMinimizer *min = new RooMinimizer( *nll_sim );

  // call the minimiser and save the result
  min->hesse();
  min->migrad();
  RooFitResult *result = min->save();

  TIterator *floatit = result->floatParsFinal().createIterator();
  TIterator *constit = result->constPars().createIterator();
  RooRealVar *floatpar, *constpar;

  while ( ( floatpar = (RooRealVar*)floatit->Next() ) ) {
    RooRealVar *init = (RooRealVar*)pars->find( floatpar->GetName() );
    cout << "NOTICE contour: fitted parameter: " << floatpar->GetName() << ", true and fitted: " 
	 << init->getVal() << "\t" << floatpar->getVal() << endl;
  }

  while ( ( constpar = (RooRealVar*)constit->Next() ) ) {
    cout << "NOTICE contour: fixed parameter: " << constpar->GetName() << endl;
  }  

  Double_t cl90 = TMath::Sqrt( ROOT::Math::chisquared_quantile(0.9,2) );

  RooPlot* contplot = min->contour( *(fu->GetVar("SinsqTh23")), *(fu->GetVar("Dm31")), 0, cl90);

  //=====================================================================================
  // set float parameters to the best-fit value and draw LLH scans
  //=====================================================================================

  while ( ( floatpar = (RooRealVar*)floatit->Next() ) ) {
    fu->GetVar( floatpar->GetName() )->setVal( floatpar->getVal() );
  }  
  
  TLegend leg1(0.3, 0.6, 0.7, 0.9);

  // optinally plot all LLH curves
  if (llhScans) {

    Int_t i = 0;
    for (auto N: nlls) {
      N.second->plotOn( &frame, LineColor(1+i), ShiftToZero(), LineStyle(1+i), Name( N.first ) );
      cout << "================================================================================" << endl;
      cout << "NOTICE contour: finished scanning channel " << N.first << endl;
      cout << "================================================================================" << endl;
      i++;
    }

    // add to legend
    for (auto N: nlls) {
      leg1.AddEntry( frame.findObject( N.first ), N.first, "l" );
    }
    leg1.SetLineWidth(0);
    leg1.SetFillStyle(0);

  }

  // clone the graphs from RooPlot for easier draw option manipulation
  TGraph *g_cont = (TGraph*)contplot->getObject(1)->Clone("g_cont");

  g_cont->SetLineColor(kRed);
  g_cont->SetMarkerColor(kRed);
  g_cont->SetLineWidth(2);

  g_cont->GetXaxis()->SetTitle("sin^2#theta_{23}");
  g_cont->GetYaxis()->SetTitle("#Delta m_{31}^2");
  g_cont->GetYaxis()->SetRangeUser(0.5*1e-3, 4.5*1e-3);
  g_cont->GetXaxis()->SetRangeUser(0., 1.);

  TLegend leg2(0.6, 0.6, 0.9, 0.9);
  leg2.AddEntry(g_cont, "2#sigma contour, sim fit", "l");
  leg2.SetLineWidth(0);
  leg2.SetFillStyle(0);

  //------------------------------------------------------------
  // write the plots to output for further manipulation
  //------------------------------------------------------------

  TFile fout((TString)outfile, "RECREATE");
  g_cont->Write("contour_90cl");
  leg2   .Write("legend_contour_90cl");
  frame  .Write("llhscan");
  leg1   .Write("legend_llhscan");
  pars  ->Write("data_parameters");  
  result->Write("fitresult");
  for (auto e: exps) e.second->Write((TString)e.first);
  fout.Close();

  //------------------------------------------------------------
  // clean-up
  //------------------------------------------------------------
  for (auto d: rfdata) if (d) delete d;
  if (min) delete min;
    
  cout << "Total run time: " << (Double_t)timer.RealTime() << endl;

}
