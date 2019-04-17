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
  // parameter randomisation
  //=====================================================================================

  vector<TString> oscpars = {"SinsqTh12", "SinsqTh13", "SinsqTh23", "dcp", "Dm21", "Dm31"};

  if (RandomisePars) {

    // set to NuFit 4.0 values and limits
    if (InvertedOrdering) o7.Set_NuFit_4p0_IO(fu);
    else o7.Set_NuFit_4p0_NO(fu);

    RooArgSet pars = fu->GetSet();
    RooArgSet obs  = fu->GetObs();
    TRandom  *rand = RooRandom::randomGenerator();
    rand->SetSeed(seed);

    TIterator *it = pars.createIterator();
    RooRealVar* var;
    while ( ( var = (RooRealVar*)it->Next() ) ) {
      
      if ( obs.find(var->GetName()) != NULL ) continue; // ignore observables

      Bool_t isOscPar = kFALSE;
      for (auto o: oscpars) isOscPar = isOscPar || ( o == (TString)var->GetName() );

      if ( isOscPar ) {
	Double_t mean  = var->getVal();
	Double_t sigma = ( var->getMax() - var->getMin() )/3 ;
	var->setVal( rand->Gaus( mean, sigma ) );
      }
      else {

	Double_t mean  = var->getVal();
	Double_t sigma = 0.2;
	var->setVal( rand->Gaus(mean, sigma) );

      }

    }

  }

  //=====================================================================================
  // parameter randomisation
  //=====================================================================================

  // create expectation value data
  std::map< string, TH1* > exps; // map with response name <--> expectation value histogram

  for (auto P: pdfs) {
    TString hname = "exp_" + P.first;
    FitPDF *pdf = P.second;
    TH3D   *exp = pdf->GetExpValHist();
    exp->SetNameTitle( hname, hname );
    exps.insert( std::make_pair( (string)P.first, (TH1*)exp ) );
  }

  // save the parameter values
  RooArgSet *pars = (RooArgSet*)fu->GetSet().snapshot(kTRUE);
  
  //=====================================================================================
  // set oscillation parameters back to central values, fix osc pars except for tm23, dm31; set
  // the systematics back to central values
  //=====================================================================================
  fu->SetNOcentvals();
  fu->GetVar("SinsqTh12")->setConstant(kTRUE);
  fu->GetVar("SinsqTh13")->setConstant(kTRUE);
  fu->GetVar("dcp")->setConstant(kTRUE);
  fu->GetVar("Dm21")->setConstant(kTRUE);

  // if systematics are not fixed set them to starting values for the fit
  if ( !fixSystematics ) {

    fu->GetVar("E_tilt")     ->setVal( 0. );
    fu->GetVar("ct_tilt")    ->setVal( 0. );
    fu->GetVar("skew_mu_amu")->setVal( 0. );
    fu->GetVar("skew_e_ae")  ->setVal( 0. );
    fu->GetVar("skew_mu_e")  ->setVal( 0. );
    fu->GetVar("NC_norm")    ->setVal( 1. );
    fu->GetVar("Tau_norm")   ->setVal( 1. );
    fu->GetVar("E_scale")    ->setVal( 0. );

  }

  // debugging
  fu->GetVar("Tau_norm")->setConstant( kTRUE );
  fu->GetVar("E_scale") ->setConstant( kTRUE );

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

  RooPlot* frame = fu->GetVar("SinsqTh23")->frame( Range(0.25, 0.75), Title("-log(L) scan vs sinsqth23") );

  vector< RooDataHist* > rfdata;         // vector to store pointers to RooFit data, for clean-up
  std::map< TString, RooAbsReal* > nlls; // map with response name <--> nll variable
  for (auto &P: pdfs) {

    TH3D   *hin = (TH3D*)exps[ (string)P.first ];
    FitPDF *pdf = P.second;
    
    TString hname = TString("rf_") + (TString)hin->GetName();
    RooDataHist *rfh = new RooDataHist( hname, hname, fu->GetObs(), Import(*hin) );

    RooAbsReal *nll = pdf->createNLL( *rfh, ExternalConstraints( o7.fPriors ), NumCPU(ncpu) );

    nlls.insert( std::make_pair( P.first, nll )  );
    rfdata.push_back( rfh );
  }

  // create a LLH scanner from simultaneous fitting, add to map
  RooAbsReal *nll_sim = simPdf.createNLL( comb, ExternalConstraints( o7.fPriors ), NumCPU(ncpu) );
  nlls.insert( std::make_pair("sim", nll_sim) );

  //=====================================================================================
  // Minimize and create contour plot
  //=====================================================================================

  if ( nlls.find(channel) == nlls.end() ) {
    cout << "Supported channels: " << endl;
    for (auto kv: nlls) cout << kv.first << endl;
    throw std::invalid_argument("ERROR! contours_wsyst unknown channel " + (string)channel);
  }

  RooMinimizer *min = new RooMinimizer( *nlls[channel] );

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
    N.second->plotOn( frame, LineColor(1+i), ShiftToZero(), LineStyle(1+i), Name( N.first ) );
    cout << "================================================================================" << endl;
    cout << "NOTICE contours_wsyst finished scanning channel " << N.first << endl;
    cout << "================================================================================" << endl;
    i++;
  }

  // add a legend
  TLegend *leg1 = new TLegend(0.3, 0.6, 0.7, 0.9);
  for (auto N: nlls) {
    leg1->AddEntry( frame->findObject( N.first ), N.first, "l" );
  }
  leg1->SetLineWidth(0);
  leg1->SetFillStyle(0);

  // draw
  TCanvas *c1 = new TCanvas("c1","c1",1);
  frame->Draw();
  leg1->Draw();

  // clone the graphs from RooPlot for easier draw option manipulation
  TGraph *g_cont = (TGraph*)contplot->getObject(1)->Clone("g_cont");

  g_cont->SetLineColor(kRed);
  g_cont->SetMarkerColor(kRed);
  g_cont->SetLineWidth(2);

  g_cont->GetXaxis()->SetTitle("sin^2#theta_{23}");
  g_cont->GetYaxis()->SetTitle("#Delta m_{31}^2");
  g_cont->GetYaxis()->SetRangeUser(0.5*1e-3, 4.5*1e-3);

  TLegend *leg2 = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg2->AddEntry(g_cont, "2#sigma contour, sim fit", "l");
  leg2->SetLineWidth(0);
  leg2->SetFillStyle(0);

  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->SetTicks();
  g_cont->Draw("AL");
  leg2->Draw();

  //------------------------------------------------------------
  // write the plots to output for further manipulation
  //------------------------------------------------------------
  TFile fout(outname,"RECREATE");
  g_cont->Write("contour_2sigma");
  leg2  ->Write("legend_contour_2sigma");
  frame ->Write("llhscan");
  leg1  ->Write("legend_llhscan");
  pars  ->Write("data_parameters");  
  result->Write("fitresult");
  for (auto e: exps) e.second->Write((TString)e.first);
  fout.Close();

  cout << "Total run time: " << (Double_t)timer.RealTime() << endl;

  // remove junk from compilation
  system("rm contours_wsyst_C*");

}

int main() {

  contours_wsyst();

}
