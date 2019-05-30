
//nmh headers
#include "AbsResponse.h"
#include "DetResponse.h"
#include "EvtResponse.h"
#include "FitUtilWsyst.h"
#include "SummaryEvent.h"
#include "SummaryParser.h"
#include "EffMass.h"
#include "FitPDF.h"
#include "NMHUtils.h"

//root and roofit headers
#include "RooArgList.h"
#include "TIterator.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooSimultaneous.h"
#include "RooPlot.h"
#include "RooNLLVar.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "TGraph.h"
#include "TF2.h"

//jpp headers
#include "JTools/JRange.hh"
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

//cpp headers
#include <iostream>
#include <stdexcept>

/** Namespace for functions and members */
namespace LLHSCAN {
  
  // structure to init PID bin settings
  struct PIDBIN {

    Double_t pid_min;
    Double_t pid_max;
    Double_t muon_cut;
    Double_t noise_cut;
    AbsResponse::reco reco;
    TString  name;

    PIDBIN() : pid_min(0), pid_max(0), muon_cut(0), noise_cut(0), reco(AbsResponse::shower), name("") {}
    PIDBIN(Double_t _pid_min, Double_t _pid_max, Double_t _muon_cut, Double_t _noise_cut, 
	   AbsResponse::reco _reco, TString _name) {
      pid_min   = _pid_min;
      pid_max   = _pid_max;
      muon_cut  = _muon_cut;
      noise_cut = _noise_cut;
      reco      = _reco;
      name      = _name;
    }

  };

  // logic to get all the parameter names in FitUtilWsyst
  vector<TString> GetFitParameters() {

    DetResponse dummyresp(DetResponse::track, "dummy");
    FitUtilWsyst dummyfu(1, dummyresp.GetHist3D(), 5, 50, -1, 0, 0, 1, EffMass::DUMMYFILE );

    RooArgList pars( dummyfu.GetSet() );
    TIterator *it = pars.createIterator();
    RooRealVar *var;

    vector< TString > parameters;
    while ( ( var = (RooRealVar*)it->Next() ) ) {
      if ( dummyfu.GetObs().find( var->GetName() ) ) continue; // ignore observables
      parameters.push_back( (TString)var->GetName() );
    }
   
    return parameters;
  }

  // functions to  be able to use the shower energy and track direction for high-purity tracks
  Double_t CustomEnergy(SummaryEvent* evt) {

    if ( evt->Get_shower_ql0() > 0.5 ) { return evt->Get_shower_energy(); }
    else                               { return evt->Get_track_energy();  }

  }

  TVector3 CustomDir(SummaryEvent *evt) { return evt->Get_track_dir();      }
  TVector3 CustomPos(SummaryEvent *evt) { return evt->Get_track_pos();      }
  Double_t CustomBY (SummaryEvent *evt) { return evt->Get_track_bjorkeny(); }

};

using namespace LLHSCAN;
using namespace RooFit;

/**
   This routine splits the data to 3 PID bins, creates expectation value data at NO central values and creates a likelihood function with respect to the expectation value data. It scans the likelihood with respect to the expectation value data in the variable specified by the user while keeping other parameters fixed or free.
*/

int main(const int argc, const char **argv) {

  TStopwatch timer;

  //======================================================
  // set some variables for command line parsing
  //======================================================

  TString  pseudo_data = "pseudodata"; // variable to identify the use of pseudo-data
  Double_t def_range   = 999.;         // variable to identify the use of default range
  vector<TString> parameters = GetFitParameters();

  vector<double> _pid_cuts   = {0., 0.3, 0.7, 1.0};  // default PID cut locations
  vector<double> _muon_cuts  = {0.03, 0.03, 0.03};   // default atm. muon suppression cuts
  vector<double> _noise_cuts = {0.03, 0.03, 0.03};   // default noise suppression cuts
  vector<string> _reco       = {"shw","shw","comb"}; // default reco's

  std::map<string, AbsResponse::reco> recomap = { {"mc"  , AbsResponse::mc_truth},
						  {"trk" , AbsResponse::track},
						  {"shw" , AbsResponse::shower},
						  {"comb", AbsResponse::customreco} };

  //======================================================
  // parse command line arguments
  //======================================================  

  TString  output_name;
  TString  meff_file;
  TString  data_file;
  JTOOLS::JRange<Double_t> par_range;
  Int_t    npoints;
  TString  par_name;
  vector<string> par_names;
  vector<Double_t> par_values;
  Int_t    ncpu;
  Double_t run_time;
  vector<double> pid_cuts;
  vector<double> muon_cuts;
  vector<double> noise_cuts;
  vector<string> reco;
  vector<string> fixed_pars;
  Bool_t profileLLH;
  Bool_t IO;
  Bool_t evtResp;

  try {

    JParser<> zap("This routine creates a likelihood scan in the specified variable to illustrate the fit sensitivity to the variable");

    zap['c'] = make_field(ncpu, "Number of CPUs used for LLH calculation") = 1;

    zap['o'] = make_field(output_name, "File where output plots are written") = "llh_scanner_out.root";
    zap['m'] = make_field(meff_file, "ROOT file with effective mass histograms, created with `apps/effective_mass` applications that use the class `common_software/EffMass`") = EffMass::DUMMYFILE;
    zap['d'] = make_field(data_file, "Summary data file in MONA format") = pseudo_data;

    zap['n'] = make_field(npoints, "Number of points in the range the LLH is evaluated. If profile LLH is performed, this becomes equivalent to the number of fits") = 30;
    zap['i'] = make_field(IO, "Inverted ordering");
    zap['t'] = make_field(run_time, "Operation time of the detector in years") = 3;
    zap['p'] = make_field(par_name, "Parameter name that is LLH scanned") = parameters;
    zap['r'] = make_field(par_range, "Parameter range for LLH scan") = JTOOLS::JRange<Double_t>(-def_range, def_range);
    
    zap['x'] = make_field(par_names, "Ordered names of the parameters for which values are provided trough argument '-y'") = std::vector<string> {};
    zap['y'] = make_field(par_values, "Ordered parameter values for parameter names specified through argument '-x'") = std::vector<Double_t> {};

    zap['P'] = make_field(pid_cuts, "PID cut locations, e.g. -P 0.+0.6+1.0 for 2 bins [0, 0.6), [0.6, 1]. By default three bins are used.") = std::vector<double>{};
    zap['M'] = make_field(muon_cuts, "Cuts for atm. muon suppression for each PID bin, by default 0.03 for default PID bins") = std::vector<double>{};
    zap['N'] = make_field(noise_cuts, "Cuts for muon suppression for each PID bin, by default 0.03 for default PID bins.") = std::vector<double>{};
    zap['R'] = make_field(reco, "Reco type (mc, shw, trk or comb) for each PID bin, 'comb' means track direction and shower energy. By default shw shw comb for default PID bins.") = std::vector<string>{};

    zap['w'] = make_field(profileLLH, "In addition to LLH scans, create a profile LLH scan. This will take more time, as at each point in the plot a fit is performed");
    zap['z'] = make_field(fixed_pars, "Names of the parameters to be fixed - this only affects profile LLH scans") = std::vector<string> {};
    
    zap['e'] = make_field(evtResp, "Use event-by-event detector response (EvtResponse), instead of binned response (DetResponse). Note that this is approx 10 times slower. EvtResponse does not work with pseudodata, so a datafileneeds to be specified. The use of EvtResponse does not require an effective mass file, this can remain default.");

    if ( zap.read(argc, argv) != 0 ) return 1;
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  if (evtResp) { cout << "NOTICE LLH_scanner: using event-by-event detector response" << endl; }
  else         { cout << "NOTICE LLH_scanner: using binned detector response" << endl; }

  // EvtResponse does not work with pseudo-data
  if (evtResp && data_file == pseudo_data) {
    throw std::invalid_argument("ERROR! LLH_scanner: EvtResponse requires actual ORCA MC summary data");
  }

  // Check that datafile and effective mass file match
  if ( data_file != pseudo_data && meff_file != EffMass::DUMMYFILE ) {
    if ( !NMHUtils::DatatagMatch( data_file, meff_file) ) throw std::invalid_argument("ERROR! LLH_scanner: datatag mismatch between the input MC summary file and the effective mass file.");
  }

  // if not specified use default values
  if ( pid_cuts.size()   == 0 ) pid_cuts   = _pid_cuts; 
  if ( muon_cuts.size()  == 0 ) muon_cuts  = _muon_cuts; 
  if ( noise_cuts.size() == 0 ) noise_cuts = _noise_cuts; 

  if ( par_names.size() != par_values.size() ) {
    throw std::invalid_argument("ERROR! LLH_scanner: specified nr of parameter names " + to_string( par_names.size() ) + " does not equal the specified nr of parameter values " + to_string( par_values.size() ) );
  }

  for (auto p: par_names) {
    TString pname = p;
    if ( std::find( parameters.begin(), parameters.end(), pname ) == parameters.end() ) {
      throw std::invalid_argument("ERROR! LLH_scanner: trying to specify (-x) a value for an unkown parameter " + p);
    }
  }

  for (auto p: fixed_pars) {
    TString pname = p;
    if ( std::find( parameters.begin(), parameters.end(), pname ) == parameters.end() ) {
      throw std::invalid_argument("ERROR! LLH_scanner: trying to fix (-z) an unkown parameter " + p);
    }
  }
  
  if ( reco.size() == 0) reco = _reco;
  else {

    for (auto r: reco) {
      if ( recomap.find(r) == recomap.end() ) {
	cout << "Known reconstruction types: " << endl;
	for (auto kv: recomap) cout << kv.first << endl;
	throw std::invalid_argument("ERROR! LLH_scanner: unknown reco type " + r);
      }

    }

  }

  if ( ( muon_cuts.size() != noise_cuts.size() ) || ( (pid_cuts.size()-1) != muon_cuts.size() ) ||
       ( muon_cuts.size() != reco.size() ) ) {

    cout << pid_cuts.size()-1 << " PID ranges:" << endl;
    for (UInt_t i = 0; i < pid_cuts.size()-1; i++) {
      cout << pid_cuts[i] << "-" << pid_cuts[i+1] << "\t";
    }
    cout << endl;
    
    cout << muon_cuts.size() << " muon cuts:" << endl;
    for (auto m: muon_cuts) cout << m << "\t";
    cout << endl;

    cout << noise_cuts.size() << " noise cuts:" << endl;
    for (auto n: noise_cuts) cout << n << "\t";
    cout << endl;

    cout << reco.size() << " reco's:" << endl;
    for (auto r: reco) cout << r << "\t";
    cout << endl;
			      
    throw std::invalid_argument("ERROR! LLH_scanner: wrong configuration of PID bins and associated muon and noise cuts. For example, to define 3 PID bins [0, 0.3), [0.3, 0.7), [0.7, 1) with corresponding muon and noise cuts, do -P 0. -P 0.3 -P 0.7 -P 1. -M 0.02 -M 0.03 -M 0.02 -N 0.01 -N 0.01 -N 0.01 -R shw -R shw -R trk");

  }

  // read the data to PIDBIN structure for further manipulation
  vector< PIDBIN > pidbins;

  for (UInt_t i = 0; i < pid_cuts.size()-1; i++) {

    if ( pid_cuts[i] >= pid_cuts[i+1] ) {
      throw std::invalid_argument("ERROR! PID cuts need to be given in increasing order");
    }
    
    TString name = "pid_" + to_string(pid_cuts[i]).substr(0,3) + "_" + to_string(pid_cuts[i+1]).substr(0,3);

    cout << "==========================================================================="<< endl;
    cout << "PID range: " << pid_cuts[i] << "-" << pid_cuts[i+1] << endl;
    cout << "muon cut : " << muon_cuts[i] << endl;
    cout << "noise cut: " << noise_cuts[i] << endl;
    cout << "reco     : " << reco[i] << endl;
    cout << "name     : " << name << endl;
    
    pidbins.push_back(PIDBIN(pid_cuts[i], pid_cuts[i+1], muon_cuts[i], noise_cuts[i], recomap[ reco[i] ], name));
  }
  cout << "==========================================================================="<< endl;
  
  //======================================================
  // initialise responses and pdf's, create data
  //======================================================

  Int_t    nebins   = 24;
  Int_t    nctbins  = 40;
  Int_t    nbybins  =  1;

  // init the responses and fill
  //---------------------

  std::map< TString, AbsResponse* > resps;

  for (UInt_t i = 0; i < pidbins.size(); i++) {

    PIDBIN PB = pidbins[i];
    TString name = PB.name+"R";
    AbsResponse *resp;

    if ( evtResp ) {
      resp = new EvtResponse( PB.reco, PB.name+"R", nebins, 1, 100, nctbins, -1, 1, nbybins, 0, 1 );
    }
    else {
      resp = new DetResponse( PB.reco, PB.name+"R", nebins, 1, 100, nctbins, -1, 1, nbybins, 0, 1 );
    }

    // decide which reco quality must be checked
    if ( PB.reco == AbsResponse::track ) {
      resp->AddCut( &SummaryEvent::Get_track_ql0 , std::greater<double>(), 0.5, true);
    }
    else if (PB.reco == AbsResponse::customreco ) {
      resp->AddCut( &SummaryEvent::Get_track_ql0 , std::greater<double>(), 0.5, true);
      resp->SetObsFuncPtrs( &CustomEnergy, &CustomDir, &CustomPos, &CustomBY );
    }
    else if ( PB.reco == AbsResponse::shower ) {
      resp->AddCut( &SummaryEvent::Get_shower_ql0, std::greater<double>(), 0.5, true);
    }
    
    // for the last PID bin the upper cut is inclusive
    if ( i == (pidbins.size()-1) ) {
      resp->AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), PB.pid_max, true);
    }
    else {
      resp->AddCut( &SummaryEvent::Get_RDF_track_score, std::less<double>(), PB.pid_max, true);
    }

    resp->AddCut(&SummaryEvent::Get_RDF_track_score, std::greater_equal<double>(), PB.pid_min  , true);
    resp->AddCut(&SummaryEvent::Get_RDF_noise_score, std::less_equal<double>()   , PB.noise_cut, true);
    resp->AddCut(&SummaryEvent::Get_RDF_muon_score , std::less_equal<double>()   , PB.muon_cut , true);

    resps.insert( std::make_pair(PB.name, resp) );
      
  }

  if ( data_file == pseudo_data ) {
    SummaryEvent evt;
    for (Int_t i = 0; i < 5e6; i++) {
      evt.FillPseudoData();
      for (auto &R: resps) R.second->Fill( &evt );
    }
  }
  else {
    SummaryParser sp( data_file );
    for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
      if ( i % 500000 == 0 ) {
	cout << "NOTICE LLH_scanner: filling responses, event " << i << " out of " 
	     << sp.GetTree()->GetEntries() << endl;
      }
      for (auto &R: resps) R.second->Fill( sp.GetEvt(i) );
    }
  }

  cout << "NOTICE LLH_scanner: responses ready" << endl;

  // init fitutil and pdf's, create data
  //---------------------
  FitUtilWsyst fu(run_time, resps.begin()->second->GetHist3DTrue(), resps.begin()->second->GetHist3DReco(), 
		  3, 65, -1, -1e-3, 0, 1, meff_file);

  std::map< TString, FitPDF* > pdfs;
  for (auto &R: resps) {
    TString pdfname = "pdf_" + R.second->GetRespName();
    FitPDF *pdf = new FitPDF(pdfname, pdfname, &fu, R.second);
    pdfs.insert( std::make_pair( R.first, pdf ) );
  }

  Double_t scan_min, scan_max;  // extract the scanner parameter minimum and maximum

  if (IO) {
    fu.FreeParLims();
    fu.SetIOcentvals();
    fu.SetIOlims();
    scan_min = fu.GetVar(par_name)->getMin();
    scan_max = fu.GetVar(par_name)->getMax();
  }
  else {
    fu.FreeParLims();
    fu.SetNOcentvals();
    fu.SetNOlims();
    scan_min = fu.GetVar(par_name)->getMin();
    scan_max = fu.GetVar(par_name)->getMax();
  }

  fu.FreeParLims(); // free the limits
  
  // perform parameter modifications requested through the command line
  for (UInt_t i = 0; i < par_names.size(); i++) {
    fu.GetVar( par_names[i] )->setVal( par_values[i] );
  }

  // fix the requested parameters
  cout << "=====================================================================" << endl;
  for (auto fp: fixed_pars) {
    fu.GetVar( fp )->setConstant(kTRUE);
    cout << "NOTICE LLH_scanner: fixed parameter " << fp << " to value " << fu.GetVar(fp)->getVal() << endl;
  }

  // notify
  cout << "=====================================================================" << endl;
  RooArgSet set = fu.GetSet();
  TIterator *it = set.createIterator();
  RooRealVar *V;
  while ( ( V = (RooRealVar*)it->Next() ) ) {
    cout << "NOTICE LLH_scanner: data created at parameter " << V->GetName() << " value " << V->getVal() << endl;
  }
  cout << "=====================================================================" << endl;
  
  RooArgSet *pars_data = (RooArgSet*)fu.GetSet().snapshot(kTRUE);

  std::map< string, TH1* > exphists;
  std::map< string, TH1* > exphists_err;
  for (auto &P: pdfs) {
    TH3D* h    = P.second->GetExpValHist();
    TH3D* herr = P.second->GetExpValErrHist();
    exphists.insert( std::make_pair( (string)P.first, (TH1*)h ) );
    exphists_err.insert( std::make_pair( (string)P.first, (TH1*)herr ) );
  }

  RooCategory categs("categs", "data categories for each data pdf pair");
  for (auto &E: exphists) {
    categs.defineType( (TString)E.first );
  }

  RooDataHist combData("combData", "combData", fu.GetObs(), categs, exphists);
  
  RooSimultaneous simPdf("simPdf","simPdf", categs);
  for (auto &P: pdfs) {
    simPdf.addPdf( *P.second, P.first );
  }

  //======================================================
  // manipulate the scanned parameter
  //======================================================

  // if limits not specified on command-line, use the limits as defined in the fit utility for IO or NO
  RooRealVar *var = fu.GetVar(par_name);
  if ( par_range.getLowerLimit() == -def_range ) par_range.setLowerLimit( scan_min );
  if ( par_range.getUpperLimit() ==  def_range ) par_range.setUpperLimit( scan_max );

  // set the parameter limits
  var->setMin( par_range.getLowerLimit() );
  var->setMax( par_range.getUpperLimit() );

  //======================================================
  // create the LLH scanners, for each PID bin and simultaneous
  //======================================================
  
  vector<RooDataHist*> rfdata; // pointers for clean-up
  std::map<TString, RooNLLVar*> nlls;

  // create NLL scanner in each PID region
  for (auto &P: pdfs) {

    TH3D   *h   = (TH3D*)exphists[ (string)P.first ];
    FitPDF *pdf = P.second;

    TString rfhname = "rf_"+P.first;
    RooDataHist *rfh = new RooDataHist(rfhname, rfhname, fu.GetObs(), Import( *h ) );

    TString nllname = "nll_" + P.first;
    RooNLLVar *nll = new RooNLLVar( nllname, nllname, *pdf, *rfh, NumCPU(ncpu) );

    nlls.insert( std::make_pair( P.first, nll ) );
    rfdata.push_back( rfh );

  }

  // create the combined scan
  RooNLLVar *nll_comb = new RooNLLVar("nll_comb","nll_comb", simPdf, combData, NumCPU(ncpu));
  nlls.insert( std::make_pair("comb", nll_comb) );

  cout << "NOTICE LLH_scanner: LLH calculators ready" << endl;

  //======================================================
  // do the plotting, write plots to output
  //======================================================

  TString title = TString("-Log(L) scan in ") + (TString)var->GetName();

  cout << "=====================================================================" << endl;
  cout << "NOTICE LLH_scanner: scanning LLH in variable " << var->GetName() << " in range " 
       << var->getMin() << " - " << var->getMax() << ", " << " nr of points " << npoints << endl;
  cout << "=====================================================================" << endl;

  std::map<TString, TGraph*> graphs;

  // plot all LLH curves. In principle, this could be achieved with nll->plotOn(frame), but identical
  // behaviour can be obtained, much faster, with just doing nll->getVal(). There seems to be some
  // integration overhead in plotOn.

  Int_t style = 0;                       // counter for line color and style
  TF2 shifter("shifter", "y - [0]", 1);  // function to help shift the Y values to start from 0

  for (auto N: nlls) {

    // create graph
    TString name = "llhscan_" + N.first;
    TGraph *gnll = new TGraph();
    gnll->SetNameTitle( name, title + ", " + N.first );

    // save the default value
    Double_t X_default = var->getVal();

    // go over scan range
    for (Double_t X = var->getMin(); X <= var->getMax(); X += ( var->getMax() - var->getMin() )/npoints ) {
      var->setVal(X);
      gnll->SetPoint( gnll->GetN(), X, N.second->getVal() );
    }

    // find the minimum Y value and shift the graph to start from 0
    shifter.SetParameter(0, TMath::MinElement( gnll->GetN(), gnll->GetY() ) );
    gnll->Apply( &shifter );

    // set the parameter back to default
    var->setVal( X_default );

    // set the style
    gnll->SetLineColor( 1+style );
    gnll->SetLineStyle( 1+style );
    gnll->SetLineWidth(2);

    // store the graph
    graphs.insert( std::make_pair(N.first, gnll) );

    // increment the style
    style++;

    cout << "=====================================================================" << endl;
    cout << "NOTICE LLH_scanner: scan for PID range " << N.first << " done" << endl;
    cout << "=====================================================================" << endl;

  } // end loop over NLL functions

  // add a legend
  TLegend *leg1 = new TLegend(0.3, 0.6, 0.7, 0.9);
  for (auto N: nlls) {
     leg1->AddEntry( graphs[N.first], N.first, "l" );
  }
  leg1->SetLineWidth(0);
  leg1->SetFillStyle(0);

  //======================================================
  // if requested, perform profile LLH scan
  //======================================================

  RooPlot frame( *var, var->getMin(), var->getMax(), npoints );
  frame.SetNameTitle(title, title);
  
  if (profileLLH) {

    // add some priors for the systematics
    RooGaussian p_mu_amu("mu_amu_prior", "mu_amu_prior", *fu.GetVar("skew_mu_amu"), RooConst(0.), RooConst(0.2) );
    RooGaussian p_e_ae  ("e_ae_prior"  , "e_ae_prior"  , *fu.GetVar("skew_e_ae")  , RooConst(0.), RooConst(0.2) );
    RooGaussian p_mu_e  ("mu_e_prior"  , "mu_e_prior"  , *fu.GetVar("skew_mu_e")  , RooConst(0.), RooConst(0.2) );
    RooGaussian pescale ("escale_prior", "escale_prior", *fu.GetVar("E_scale")    , RooConst(0.), RooConst(0.15));
    RooGaussian pncnorm ("ncnorm_prior", "ncnorm_prior", *fu.GetVar("NC_norm")    , RooConst(1.), RooConst(0.1) );

    RooArgSet priors ( p_mu_amu, p_e_ae, p_mu_e, pescale, pncnorm );

    // create a new LLH function that includes the priors
    RooNLLVar *nll_comb_wp = new RooNLLVar("nll_comb_wp","nll_comb_wp", simPdf, combData, ExternalConstraints(priors), NumCPU(ncpu) );
    
    cout << "=====================================================================" << endl;
    cout << "NOTICE LLH_scanner: started profile likelihood calculation" << endl;
    cout << "=====================================================================" << endl;
    RooAbsReal *profile = nll_comb_wp->createProfile( *var );
    profile->plotOn( &frame, LineColor(style+1), LineStyle(style+1), ShiftToZero(), Name("profileLLH") );
    leg1->AddEntry( frame.findObject( "profileLLH" ), "profileLLH", "l" );
     
  }
  
  TFile fout(output_name, "RECREATE");
  frame.Write("LLH_plot");
  leg1->Write("LLH_legend");
  pars_data->Write("Parameters");
  for (auto E: exphists)     E.second->Write( (TString)E.first );
  for (auto E: exphists_err) E.second->Write( (TString)E.first + TString("_err") );
  for (auto g: graphs)       g.second->Write();
  fout.Close();

  // clean-up
  for (auto R: resps) if (R.second) delete R.second;
  for (auto P: pdfs) if (P.second) delete P.second;
  for (auto E: exphists) if (E.second) delete E.second;
  for (auto E: exphists_err) if (E.second) delete E.second;
  for (auto d: rfdata) if (d) delete d;
  for (auto n: nlls) if (n.second) delete n.second;
  for (auto g: graphs) if (g.second) delete g.second;
  
  cout << "NOTICE LLH_scanner: total run time [s]: " << (Double_t)timer.RealTime() << endl;

}
