#include "Asymmetry.h"
#include "NMHUtils.h"
#include "FileHeader.h"

#include "TFile.h"
#include "TMath.h"

#include <iostream>
#include <stdexcept>

using namespace std;
using namespace ASYM;

/**
   Asymmetry analysis.

   \param summary_files      List of all available summary files
   \parm  NH_IH_pairs_list   List with NH and IH file pair on each line (NH first, IH second), separated by a space
   \param idstr              ID string added to the output files
 */
void Asymmetry(TString summary_files, TString NH_IH_pairs_list, TString idstr="") {


  // get the flux histograms from the FluxChain.C output files
  vector<TString> pair_lines = NMHUtils::ReadLines(NH_IH_pairs_list);
  vector< std::pair<TString, TString> > NH_IH_pairs;

  for (auto &l: pair_lines) {
    TString first = l( 0, l.Index(" ") );
    TString second = l( l.Index(" "), l.Length() );
    NH_IH_pairs.push_back( std::make_pair(first, second) );
  }

  InitFluxHists(NH_IH_pairs);
  
  // initialize other variables/hists and fill the histograms for weighting
  InitVars(std::get<NH_nu>(fFluxHists[0][0]), summary_files);
  FillSimHists();

  // init event selections, fill them and calculate the asymmetries
  InitEvtSels_1( NH_IH_pairs.size() );
  FillSelections();
  GetAsymmetries(NH_IH_pairs, idstr);

  // cleanup
  CleanUp();

}

//*****************************************************************************************

void ASYM::InitFluxHists( vector< std::pair<TString, TString> > &NH_IH_pairs ) {

  // loop over NH and IH file pairs
  for (auto &p: NH_IH_pairs) {

    TFile nhfile(p.first , "READ");
    TFile ihfile(p.second, "READ");

    // per each pair there is a vector with histograms per flavor
    vector< std::tuple<TH2D*,TH2D*,TH2D*,TH2D*> > p_hists;

    // loop over flavors
    for (Int_t f = elec; f <= nc; f++) {
      
      TH2D *h_NH_det_nu, *h_NH_det_nub, *h_IH_det_nu, *h_IH_det_nub;

      // set the pointers to the correct histograms
      if (f < nc) {
	
	TString hname_nu  = "detflux/detflux_" + fFlav_str[f] + "_cc_nu";
	TString hname_nub = "detflux/detflux_" + fFlav_str[f] + "_cc_nub";

	h_NH_det_nu  = (TH2D*)nhfile.Get(hname_nu);
	h_NH_det_nub = (TH2D*)nhfile.Get(hname_nub);
	h_IH_det_nu  = (TH2D*)ihfile.Get(hname_nu);
	h_IH_det_nub = (TH2D*)ihfile.Get(hname_nub);

      }

      else {

	h_NH_det_nu = (TH2D*)nhfile.Get("detflux/detflux_elec_nc_nu");
	h_NH_det_nu->Add( (TH2D*)nhfile.Get("detflux/detflux_muon_nc_nu") );
	h_NH_det_nu->Add( (TH2D*)nhfile.Get("detflux/detflux_tau_nc_nu") );

	h_NH_det_nub = (TH2D*)nhfile.Get("detflux/detflux_elec_nc_nub");
	h_NH_det_nub->Add( (TH2D*)nhfile.Get("detflux/detflux_muon_nc_nub") );
	h_NH_det_nub->Add( (TH2D*)nhfile.Get("detflux/detflux_tau_nc_nub") );

	h_IH_det_nu = (TH2D*)ihfile.Get("detflux/detflux_elec_nc_nu");
	h_IH_det_nu->Add( (TH2D*)ihfile.Get("detflux/detflux_muon_nc_nu") );
	h_IH_det_nu->Add( (TH2D*)ihfile.Get("detflux/detflux_tau_nc_nu") );

	h_IH_det_nub = (TH2D*)ihfile.Get("detflux/detflux_elec_nc_nub");
	h_IH_det_nub->Add( (TH2D*)ihfile.Get("detflux/detflux_muon_nc_nub") );
	h_IH_det_nub->Add( (TH2D*)ihfile.Get("detflux/detflux_tau_nc_nub") );

      }

      // per each flavor there are 4 histograms
      auto f_hists = std::make_tuple( (TH2D*)h_NH_det_nu->Clone(), (TH2D*)h_NH_det_nub->Clone(),
				      (TH2D*)h_IH_det_nu->Clone(), (TH2D*)h_IH_det_nub->Clone() );

      std::get<NH_nu> (f_hists)->SetDirectory(0);
      std::get<NH_nub>(f_hists)->SetDirectory(0);
      std::get<IH_nu> (f_hists)->SetDirectory(0);
      std::get<IH_nub>(f_hists)->SetDirectory(0);

      // add the histograms to the flavors vector
      p_hists.push_back( f_hists );
      
    }

    // add the flavors vector to the vector for pairs
    fFluxHists.push_back(p_hists);

  }

  cout << "NOTICE ASYM::InitFluxHists() finished initalising flux histograms." << endl;

}

//*****************************************************************************************

void ASYM::InitVars(TH2D *h_template, TString summary_files) {

  // init summary parser and add all files for reading

  vector<TString> sfiles = NMHUtils::ReadLines(summary_files);
  fSp = new SummaryParser(sfiles[0]);
  for (UInt_t i = 1; i < sfiles.size(); i++) { fSp->GetTree()->Add( sfiles[i] ); }


  // initialize weight histogram for each flavor

  for (Int_t f = elec; f <= nc; f++) {

    TH2D* h_nu  = (TH2D*)h_template->Clone();
    TH2D* h_nub = (TH2D*)h_template->Clone();
    h_nu->Reset();
    h_nub->Reset();

    TString hname_nu  = "h_w_" + fFlav_str[f] + "_nu";
    TString hname_nub = "h_w_" + fFlav_str[f] + "_nub";
    h_nu ->SetNameTitle(hname_nu , hname_nu);
    h_nub->SetNameTitle(hname_nub, hname_nub);
    h_nu->SetDirectory(0);
    h_nub->SetDirectory(0);

    fSimHists.push_back( std::make_pair(h_nu, h_nub) );

  }

  cout << "NOTICE ASYM::InitVars() finished initalising variables." << endl; 
 
}

//*****************************************************************************************

void ASYM::FillSimHists() {

  for (Int_t i = 0; i < fSp->GetTree()->GetEntries(); i++) {
    
    fSp->GetTree()->GetEntry(i);
    SummaryEvent *evt = fSp->GetEvt();
    
    // ignore muons and noise
    if ( evt->Get_MC_is_neutrino() < 0.5 ) continue;

    // charged current
    if ( evt->Get_MC_is_CC() > 0.5 ) {

      // check for the true flavor of the neutrino
      for (Int_t f = elec; f <= tau; f++) {

	if ( TMath::Abs( evt->Get_MC_type() ) == fFlav_pdg[f] ) {

	  // fill hist, depending on whether it is nu/nub
	  if ( evt->Get_MC_type() > 0 ) {
	    fSimHists[f].first ->Fill( evt->Get_MC_energy(), -evt->Get_MC_dir_z() );
	  }
	  else {
	    fSimHists[f].second->Fill( evt->Get_MC_energy(), -evt->Get_MC_dir_z() );
	  }

	}

      } // end loop over flavors

    }
    // neutral current
    else {

      // fill nc hist, depending on whether it is nu/nub
      if ( evt->Get_MC_type() > 0 ) {
	fSimHists[nc].first ->Fill( evt->Get_MC_energy(), -evt->Get_MC_dir_z() );
      }
      else {
	fSimHists[nc].second->Fill( evt->Get_MC_energy(), -evt->Get_MC_dir_z() );
      }

    }

  }

  cout << "NOTICE ASYM::FillSimHists() finished filling histograms for weighting." << endl;

}

//*****************************************************************************************

void ASYM::InitEvtSels_1(Int_t NH_IH_pcount) {

  Double_t emin  =   2;
  Double_t emax  = 100;
  Double_t ebins =  40;

  Double_t ctmin  = -1;
  Double_t ctmax  =  0;
  Double_t ctbins = 40;

  EventSelection t(EventSelection::track, "track", NULL, ebins, emin, emax, ctbins, ctmin, ctmax);
  t.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
  t.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
  t.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
  t.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
  t.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

  EventSelection s(EventSelection::shower, "shower", NULL, ebins, emin, emax, ctbins, ctmin, ctmax);
  s.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
  s.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
  s.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(),  0.6, true );
  s.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  s.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

  //----------DEBUG----------------------------
  // remove downward going events for swim/paramNMH comparisons
  //-------------------------------------------
  t.AddCut( &SummaryEvent::Get_MC_dir_z       , std::greater<double>()   ,  0., true );
  s.AddCut( &SummaryEvent::Get_MC_dir_z       , std::greater<double>()   ,  0., true );

  // use MC truth energy and angle and exact PID (only muons)
  EventSelection t_mc(EventSelection::mc_truth, "track_mctruth", NULL, ebins, emin, emax, ctbins, ctmin, ctmax);
  t_mc.AddCut( &SummaryEvent::Get_track_ql0      , std::greater<double>()   ,  0.5, true );
  t_mc.AddCut( &SummaryEvent::Get_track_ql1      , std::greater<double>()   ,  0.5, true );
  t_mc.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  t_mc.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), 0.18, true );
  t_mc.AddCut( &SummaryEvent::Get_MC_type        , std::equal_to<double>()  ,   14, false );
  t_mc.AddCut( &SummaryEvent::Get_MC_type        , std::equal_to<double>()  ,  -14, false );

  // use MC truth energy and angle and exact PID (everything but muons)
  EventSelection s_mc(EventSelection::mc_truth, "shower_mctruth", NULL, ebins, emin, emax, ctbins, ctmin, ctmax);
  s_mc.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
  s_mc.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
  s_mc.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
  s_mc.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );
  s_mc.AddCut( &SummaryEvent::Get_MC_type        , std::not_equal_to<double>(),   14, true );
  s_mc.AddCut( &SummaryEvent::Get_MC_type        , std::not_equal_to<double>(),  -14, true );

  for (Int_t i = 0; i < NH_IH_pcount; i++) {
   
    vector<EventSelection*> NH_sels = { new EventSelection(t), new EventSelection(s), new EventSelection(t_mc), new EventSelection(s_mc) };
    vector<EventSelection*> IH_sels = { new EventSelection(t), new EventSelection(s), new EventSelection(t_mc), new EventSelection(s_mc) };
    fEvtSels.push_back( std::make_pair(NH_sels, IH_sels) );

  }

  cout << "NOTICE ASYM::InitEvtSels_1() finished initalising event selections." << endl;

}

//*****************************************************************************************

void ASYM::FillSelections() {

  for (Int_t i = 0; i < fSp->GetTree()->GetEntries(); i++) {

    fSp->GetTree()->GetEntry(i);
    SummaryEvent *evt = fSp->GetEvt();

    if ( !evt->Get_MC_is_neutrino() ) continue;
      
    // determine the flavor of this event
    Int_t flav = -1;

    if ( !evt->Get_MC_is_CC() ) {
      flav = nc;
    }
    else {

      for (Int_t f = 0; f <= tau; f++) {
	if ( TMath::Abs( evt->Get_MC_type() ) == fFlav_pdg[f] ) {
	  flav = f;
	  break;
	}
      }

    }

    // loop over NH/IH pairs
    for (UInt_t p = 0; p < fFluxHists.size(); p++) {
      
      // get the detected flux histograms for this NH/IH pair and the simulated hists
      TH2D *h_NH_det_nu  = std::get<NH_nu> (fFluxHists[p][flav]);
      TH2D *h_NH_det_nub = std::get<NH_nub>(fFluxHists[p][flav]);
      TH2D *h_IH_det_nu  = std::get<IH_nu> (fFluxHists[p][flav]);
      TH2D *h_IH_det_nub = std::get<IH_nub>(fFluxHists[p][flav]);

      TH2D* h_sim_nu  = fSimHists[flav].first;
      TH2D* h_sim_nub = fSimHists[flav].second;

      // get bins for weight calculation
      Int_t ebin  = h_NH_det_nu->GetXaxis()->FindBin(  evt->Get_MC_energy() );
      Int_t ctbin = h_NH_det_nu->GetYaxis()->FindBin( -evt->Get_MC_dir_z()  );

      // calculate weights for this event, NH and IH
      Double_t w_NH, w_IH;
      if ( evt->Get_MC_type() > 0 ) {
	w_NH = h_NH_det_nu->GetBinContent(ebin, ctbin)/h_sim_nu->GetBinContent(ebin, ctbin);
	w_IH = h_IH_det_nu->GetBinContent(ebin, ctbin)/h_sim_nu->GetBinContent(ebin, ctbin);
      }
      else {
	w_NH = h_NH_det_nub->GetBinContent(ebin, ctbin)/h_sim_nub->GetBinContent(ebin, ctbin);
	w_IH = h_IH_det_nub->GetBinContent(ebin, ctbin)/h_sim_nub->GetBinContent(ebin, ctbin);
      }

      // call fill for all event classes for NH and IH with this event
      for (auto ev: fEvtSels[p].first)  ev->Fill(evt, w_NH);
      for (auto ev: fEvtSels[p].second) ev->Fill(evt, w_IH);

    }

  }

  cout << "NOTICE ASYM::FillSelections() finished filling event selections." << endl;

}

//*****************************************************************************************

void ASYM::GetAsymmetries(vector< std::pair<TString, TString> > NH_IH_pairs, TString idstr) {

  // loop over NH/IH pairs
  for (UInt_t p = 0; p < fEvtSels.size(); p++) {

    vector<EventSelection*> NHsels = fEvtSels[p].first;
    vector<EventSelection*> IHsels = fEvtSels[p].second;

    if ( NHsels.size() != IHsels.size() ) {
      throw std::invalid_argument( " ERROR! ASYM::GetAsymmetries() unexpected event selections." );
    }

    // store the input FluxChain.C output filenames in the header
    FileHeader Asym("Asymmetry");
    Asym.AddParameter("NHflux", NH_IH_pairs[p].first);
    Asym.AddParameter("IHflux", NH_IH_pairs[p].second);
    vector<TH2D*> hout;

    // for each NH/IH pair, loop over the EventSelection classes
    for (UInt_t c = 0; c < NHsels.size(); c++) {
      
      TString selname = NHsels[c]->Get_SelName();

      auto a_data = NMHUtils::Asymmetry( NHsels[c]->Get_h_E_costh(), IHsels[c]->Get_h_E_costh(),
					 "Asymmetry_" + selname, 0, 100, -1, 0);

      // add the asymmetry, chi2 and number of deg. of freedom to the header for each class
      Asym.AddParameter("Asymmetry_" + selname, (TString)to_string( std::get<1>(a_data) ) );
      Asym.AddParameter("Chi2_" + selname     , (TString)to_string( std::get<2>(a_data) ) );
      Asym.AddParameter("Ndof_" + selname     , (TString)to_string( std::get<3>(a_data) ) );
      
      // add the asymmetry and event class histograms
      hout.push_back( std::get<0>(a_data) );
      hout.push_back( (TH2D*)NHsels[c]->Get_h_E_costh()->Clone("NH_" + selname) );
      hout.push_back( (TH2D*)IHsels[c]->Get_h_E_costh()->Clone("IH_" + selname) );

    }

    // write the output to file
    TFile fout("output/Asymmetry_sel" + (TString)to_string(p) + idstr + ".root", "RECREATE");
    Asym.WriteHeader(&fout);
    for (auto &h: hout) { h->Write(); delete h; }
    fout.Close();

  }

  cout << "NOTICE ASYM::GetAsymmetries() finished calculating asymmetries." << endl;

}

//*****************************************************************************************

void ASYM::CleanUp() {

  if (fSp) delete fSp;
  
  for (auto &p: fSimHists) {
    if (p.first) delete p.first;
    if (p.second) delete p.second;
  }

  for (auto &fv: fFluxHists) {
    for (auto &t: fv) {
      if ( std::get<NH_nu>(t) )  delete std::get<NH_nu>(t);
      if ( std::get<NH_nub>(t) ) delete std::get<NH_nub>(t);
      if ( std::get<IH_nu>(t) )  delete std::get<IH_nu>(t);
      if ( std::get<IH_nub>(t) ) delete std::get<IH_nub>(t);
    }
  }

  for (auto &p: fEvtSels) {
    for (auto &ev: p.first)  { if (ev) delete ev; }
    for (auto &ev: p.second) { if (ev) delete ev; }
  }

  cout << "NOTICE ASYM::CleanUp() finished." << endl;
}
