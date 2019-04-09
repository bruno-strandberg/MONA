
//nmh headers
#include "DetResponse.h"
#include "FitUtilWsyst.h"
#include "SummaryEvent.h"
#include "SummaryParser.h"
#include "EffMass.h"
#include "FitPDF.h"

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
    DetResponse::reco reco;
    TString  name;

    PIDBIN() : pid_min(0), pid_max(0), muon_cut(0), noise_cut(0), reco(DetResponse::shower), name("") {}
    PIDBIN(Double_t _pid_min, Double_t _pid_max, Double_t _muon_cut, Double_t _noise_cut, 
	   DetResponse::reco _reco, TString _name) {
      pid_min   = _pid_min;
      pid_max   = _pid_max;
      muon_cut  = _muon_cut;
      noise_cut = _noise_cut;
      reco      = _reco;
      name = _name;
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

  std::map<string, DetResponse::reco> recomap = { {"mc"  , DetResponse::mc_truth},
						  {"trk" , DetResponse::track},
						  {"shw" , DetResponse::shower},
						  {"comb", DetResponse::customreco} };

  //======================================================
  // parse command line arguments
  //======================================================  

  TString  output_name;
  TString  meff_file;
  TString  data_file;
  JTOOLS::JRange<Double_t> par_range;
  Int_t    npoints;
  TString  par_name;
  Int_t    ncpu;
  vector<double> pid_cuts;
  vector<double> muon_cuts;
  vector<double> noise_cuts;
  vector<string> reco;

  try {

    JParser<> zap("This routine creates a likelihood scan in the specified variable to illustrate the fit sensitivity to the variable");

    zap['c'] = make_field(ncpu, "Number of CPUs used for LLH calculation") = 1;

    zap['o'] = make_field(output_name, "File where output plots are written") = "llh_scanner_out.root";
    zap['m'] = make_field(meff_file, "ROOT file with effective mass histograms, created with `apps/effective_mass` applications that use the class `common_software/EffMass`") = EffMass::DUMMYFILE;
    zap['d'] = make_field(data_file, "Summary data file in MONA format") = pseudo_data;

    zap['p'] = make_field(par_name, "Parameter name that is LLH scanned") = parameters;
    zap['r'] = make_field(par_range, "Parameter range for LLH scan") = JTOOLS::JRange<Double_t>(-def_range, def_range);
    zap['n'] = make_field(npoints, "Number of points in the range the LLH is evaluated") = 100;

    zap['P'] = make_field(pid_cuts, "PID cut locations, e.g. -P 0.+0.6+1.0 for 2 bins [0, 0.6), [0.6, 1]. By default three bins are used.") = std::vector<double>{};
    zap['M'] = make_field(muon_cuts, "Cuts for atm. muon suppression for each PID bin, by default 0.03 for default PID bins") = std::vector<double>{};
    zap['N'] = make_field(noise_cuts, "Cuts for muon suppression for each PID bin, by default 0.03 for default PID bins.") = std::vector<double>{};
    zap['R'] = make_field(reco, "Reco type (mc, shw, trk or comb) for each PID bin, 'comb' means track direction and shower energy. By default shw shw comb for default PID bins.") = std::vector<string>{};

    if ( zap.read(argc, argv) != 0 ) return 1;
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  // if not specified use default values
  if ( pid_cuts.size()   == 0 ) pid_cuts   = _pid_cuts; 
  if ( muon_cuts.size()  == 0 ) muon_cuts  = _muon_cuts; 
  if ( noise_cuts.size() == 0 ) noise_cuts = _noise_cuts; 

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

  std::map< TString, DetResponse* > resps;

  for (UInt_t i = 0; i < pidbins.size(); i++) {

    PIDBIN PB = pidbins[i];
    TString name = PB.name+"R";
    DetResponse *resp = new DetResponse( PB.reco, PB.name+"R", nebins, 1, 100, nctbins, -1, 1, nbybins, 0, 1 );

    // decide which reco quality must be checked
    if ( PB.reco == DetResponse::track ) {
      resp->AddCut( &SummaryEvent::Get_track_ql0 , std::greater<double>(), 0.5, true);
    }
    else if (PB.reco == DetResponse::customreco ) {
      resp->AddCut( &SummaryEvent::Get_track_ql0 , std::greater<double>(), 0.5, true);
      resp->SetObsFuncPtrs( &CustomEnergy, &CustomDir, &CustomPos, &CustomBY );
    }
    else if ( PB.reco == DetResponse::shower ) {
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
  FitUtilWsyst fu(3, resps.begin()->second->GetHist3D(), 3, 65, -1, -1e-3, 0, 1, meff_file);

  std::map< TString, FitPDF* > pdfs;
  for (auto &R: resps) {
    TString pdfname = "pdf_" + R.second->GetRespName();
    FitPDF *pdf = new FitPDF(pdfname, pdfname, &fu, R.second);
    pdfs.insert( std::make_pair( R.first, pdf ) );
  }

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

  cout << "NOTICE LLH_scanner: pdf's ready" << endl;

  //======================================================
  // manipulate the scanned parameter
  //======================================================
  fu.SetNOlims();

  // if limits not specified on command-line, use the limits as defined in the fit utility
  RooRealVar *var = fu.GetVar(par_name);
  if ( par_range.getLowerLimit() == -def_range ) par_range.setLowerLimit( var->getMin() );
  if ( par_range.getUpperLimit() ==  def_range ) par_range.setUpperLimit( var->getMax() );

  // set the parameter limits
  var->setMin( par_range.getLowerLimit() );
  var->setMax( par_range.getUpperLimit() );

  //======================================================
  // create the LLH scanners, for each PID bin and simultaneous
  //======================================================
  
  TString title = TString("-Log(L) scan in ") + (TString)var->GetName();
  RooPlot frame( *var, var->getMin(), var->getMax(), npoints );
  frame.SetNameTitle(title, title);

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

  cout << "=====================================================================" << endl;
  cout << "NOTICE LLH_scanner: scanning LLH in variable " << var->GetName() << " in range " 
       << var->getMin() << " - " << var->getMax() << endl;
  cout << "=====================================================================" << endl;

  // plot all LLH curves
  Int_t i = 0;
  for (auto N: nlls) {
    N.second->plotOn( &frame, LineColor(1+i), ShiftToZero(), LineStyle(1+i), Name( N.first ) );
    cout << "=====================================================================" << endl;
    cout << "NOTICE LLH_scanner: scan for PID range " << N.first << " done" << endl;
    cout << "=====================================================================" << endl;
    i++;
  }

  // add a legend
  TLegend *leg1 = new TLegend(0.3, 0.6, 0.7, 0.9);
  for (auto N: nlls) {
     leg1->AddEntry( frame.findObject( N.first ), N.first, "l" );
  }
  leg1->SetLineWidth(0);
  leg1->SetFillStyle(0);

  TFile fout(output_name, "RECREATE");
  frame.Write("LLH_plot");
  leg1->Write("LLH_legend");
  pars_data->Write("Parameters");
  for (auto E: exphists)     E.second->Write( (TString)E.first );
  for (auto E: exphists_err) E.second->Write( (TString)E.first + TString("_err") );
  fout.Close();

  // clean-up
  for (auto R: resps) if (R.second) delete R.second;
  for (auto P: pdfs) if (P.second) delete P.second;
  for (auto E: exphists) if (E.second) delete E.second;
  for (auto E: exphists_err) if (E.second) delete E.second;
  for (auto d: rfdata) if (d) delete d;
  for (auto n: nlls) if (n.second) delete n.second;
  
  cout << "NOTICE LLH_scanner: total run time [s]: " << (Double_t)timer.RealTime() << endl;

}
