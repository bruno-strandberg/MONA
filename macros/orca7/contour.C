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
  Bool_t   RandomisePars;
  Int_t    seed;
  Int_t    ncpu;
  string   outfile;
  
  try {

    JParser<> zap("This application creates contour plots for ORCA7");

    zap['i'] = make_field(InvertedOrdering, "Inverted NMO");
    zap['t'] = make_field(sinsqth23, "SinsqTh23 value at which the 'data' is created") = 0.58;
    zap['s'] = make_field(IncludeSystematics, "Include systematic parameters");
    zap['p'] = make_field(IncludePriors, "Include priors for systematics");
    zap['r'] = make_field(RandomisePars, "Randomise all parameters when 'data' is created");
    zap['S'] = make_field(seed, "Seed for parameter randomisation") = 416;
    zap['N'] = make_field(ncpu, "Number of CPUs for minimisation") = 1;
    zap['o'] = make_field(outfile, "Output file") = "contour.root";

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
  // create vectors with only oscpars and only systpars
  //=====================================================================================

  vector<RooRealVar*> oscpars = { fu->GetVar("SinsqTh12"), fu->GetVar("SinsqTh13"), fu->GetVar("SinsqTh23"),
				  fu->GetVar("dcp"), fu->GetVar("Dm21"), fu->GetVar("Dm31") };
  vector<RooRealVar*> systpars;

  RooArgSet parset = fu->GetSet();
  TIterator *it = parset.createIterator();
  RooRealVar* var;

  while ( ( var = (RooRealVar*)it->Next() ) ) {
    
    if ( fu->GetObs().find(var->GetName()) != NULL ) continue; // ignore observables

      Bool_t isOscPar = kFALSE;
      for (auto o: oscpars) isOscPar = isOscPar || ( o == var );
      if ( !isOscPar ) systpars.push_back(var);
      
  }

  //=====================================================================================
  // randomise parameters
  //=====================================================================================

  TRandom  *rand = RooRandom::randomGenerator();
  rand->SetSeed(seed);

  // randomisation of oscillation parameters
  if ( RandomisePars ) {

    // set to NuFit 4.0 values and limits
    if (InvertedOrdering) o7.Set_NuFit_4p0_IO(fu);
    else o7.Set_NuFit_4p0_NO(fu);

    for (auto var: oscpars) {
      cout << "NOTICE contour: randomising parameter " << var->GetName() << endl; 
      Double_t mean  = var->getVal();
      Double_t sigma = ( var->getMax() - var->getMin() )/3 ;
      var->setVal( rand->Gaus( mean, sigma ) );
    }
    
  }

  // randomisation of systematics, store also the default values
  std::map<RooRealVar*, Double_t> default_systs;

  for (auto var: systpars) {

    default_systs.insert( std::make_pair( var, var->getVal() ) );
    
    if ( !IncludeSystematics ) {
      cout << "NOTICE contour: fixing parameter " << var->GetName() << endl;
      var->setConstant(kTRUE);
    }
    else  {

      if ( RandomisePars ) {
	cout << "NOTICE contour: randomising parameter " << var->GetName() << endl; 
	Double_t mean  = var->getVal();
	Double_t sigma = 0.15;
	var->setVal( rand->Gaus(mean, sigma) );

      }
    }
    
  }

  // theta-23 is set to the input value
  fu->GetVar("SinsqTh23")->setVal( sinsqth23 );

  // save the parameter values at which data is created
  RooArgSet *pars = (RooArgSet*)fu->GetSet().snapshot(kTRUE);
  
  //=====================================================================================
  // create data
  //=====================================================================================
  
  std::map< string, TH1* > exps; // map with response name <--> expectation value histogram

  for (auto P: pdfs) {
    TString hname = "exp_" + P.first;
    FitPDF *pdf = P.second;
    TH3D   *exp = pdf->GetExpValHist();
    exp->SetNameTitle( hname, hname );
    exps.insert( std::make_pair( (string)P.first, (TH1*)exp ) );
  }
  
  //=====================================================================================
  // set parameter values for fit start
  //=====================================================================================

  // for fit start, set the osc parameters to NuFit central values
  if (InvertedOrdering) { o7.Set_NuFit_4p0_IO(fu); }
  else                  { o7.Set_NuFit_4p0_NO(fu); }

  fu->FreeParLims(); // release the limits
  
  // theta-23 is set to the input value, so that the fit starts in the right quadrant
  fu->GetVar("SinsqTh23")->setVal( sinsqth23 );

  // systematics are set to default values
  for (auto var: systpars) {
    var->setVal( default_systs[var] );
  }
  
  // always fix these osc pars, no sensitivity
  fu->GetVar("SinsqTh12")->setConstant(kTRUE);
  fu->GetVar("SinsqTh13")->setConstant(kTRUE);
  fu->GetVar("dcp")->setConstant(kTRUE);
  fu->GetVar("Dm21")->setConstant(kTRUE);

  // for dev, I currently fix tau and E_scale
  // fu->GetVar("Tau_norm")->setConstant( kTRUE );
  // fu->GetVar("E_scale") ->setConstant( kTRUE );

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
    cout << "NOTICE fitted parameter: " << floatpar->GetName() << ", true and fitted: " 
	 << init->getVal() << "\t" << floatpar->getVal() << endl;
  }

  while ( ( constpar = (RooRealVar*)constit->Next() ) ) {
    cout << "NOTICE fixed parameter: " << constpar->GetName() << endl;
  }  

  Double_t cl90 = TMath::Sqrt( ROOT::Math::chisquared_quantile(0.9,2) );

  RooPlot* contplot = min->contour( *(fu->GetVar("SinsqTh23")), *(fu->GetVar("Dm31")), 0, cl90);

  //=====================================================================================
  // set float parameters to the best-fit value and draw LLH scans
  //=====================================================================================

  while ( ( floatpar = (RooRealVar*)floatit->Next() ) ) {
    fu->GetVar( floatpar->GetName() )->setVal( floatpar->getVal() );
  }  
  
  // plot all LLH curves
  Int_t i = 0;
  for (auto N: nlls) {
    N.second->plotOn( &frame, LineColor(1+i), ShiftToZero(), LineStyle(1+i), Name( N.first ) );
    cout << "================================================================================" << endl;
    cout << "NOTICE contours_wsyst finished scanning channel " << N.first << endl;
    cout << "================================================================================" << endl;
    i++;
  }

  // add a legend
  TLegend leg1(0.3, 0.6, 0.7, 0.9);
  for (auto N: nlls) {
    leg1.AddEntry( frame.findObject( N.first ), N.first, "l" );
  }
  leg1.SetLineWidth(0);
  leg1.SetFillStyle(0);

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
  g_cont->Write("contour_2sigma");
  leg2   .Write("legend_contour_2sigma");
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
