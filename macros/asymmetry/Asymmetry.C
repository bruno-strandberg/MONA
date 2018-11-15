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
   This program performs an asymmetry analysis between normal and inverted hierachy.

   The idea is that one can analyse several NH/IH pairs simultaneously, meaning that only one loop over the summary events is required.

   The program outputs files in format `output/Asymmetry_sel{pidx}_{idstr}.root`. Index `pidx` is the index of the NH/IH pair in the input file `NH_IH_pairs_list`, `idstr` is the identifier string provided to the program. Each output contains NH and IH histograms for all of the event selections in `fEvtSels`, as well as an asymmetry histogram per each event selection between NH and IH. In addition, each file contains a `FileHeader` that stores the calculated asymmetry data.

   Comment: Although the application works fine, the member variable structure is (unnecessarily?) complicated. A new asymmetry calculator based solely on the detector response should be implemented in the future.

   \param summary_files      List of all available summary files (currently noise and atm. mu are ignored), one file per line.
   \parm  NH_IH_pairs_list   List with NH and IH file pairs on each line (NH first, IH second), separated by a space Example line `../evt_sampler/output/FluxChain/file1_NH.root ../evt_sampler/output/FluxChain/file1_IH.root`. The file pairs should be created with the 'FluxChain.C' application with using effective mass info.
   \param idstr              ID string added to the output files
 */
void Asymmetry(TString summary_files, TString NH_IH_pairs_list, TString idstr="") {

  // init summary parser and add all files for reading
  vector<TString> sfiles = NMHUtils::ReadLines(summary_files);
  fSp = new SummaryParser(sfiles[0]);
  for (UInt_t i = 1; i < sfiles.size(); i++) { fSp->GetTree()->Add( sfiles[i] ); }

  // get the flux histograms from the FluxChain.C output files
  vector<TString> pair_lines = NMHUtils::ReadLines(NH_IH_pairs_list);
  vector< std::pair<TString, TString> > NH_IH_pairs;

  for (auto &l: pair_lines) {
    TString first = l( 0, l.Index(" ") );
    TString second = l( l.Index(" "), l.Length() );
    NH_IH_pairs.push_back( std::make_pair(first, second) );
  }

  // init members to be filled
  InitFluxHists(NH_IH_pairs);
  InitVars();
  InitEvtSels_1();

  // do the filling
  FillSimHists();
  FillSelections();
  FillFromResponse();

  // asymmetry calculation and writing out
  GetAsymmetries(NH_IH_pairs, idstr);

  // cleanup
  CleanUp();

}

//*****************************************************************************************

/**
   This function initialises the structure `ASYM::fFluxHists` with histograms from the `FluxChain.C` program input to the `Asymmetry.C` program.

   \param NH_IH_pairs     Vector of filename pairs, vec[index].first is the NH file, vec[index].second the IH file.
 */
void ASYM::InitFluxHists( vector< std::pair<TString, TString> > &NH_IH_pairs ) {

  // loop over NH and IH file pairs
  for (auto &p: NH_IH_pairs) {

    TFile nhfile(p.first , "READ");
    TFile ihfile(p.second, "READ");

    // per each pair there is a vector with histograms per flavor
    vector< std::tuple<TH3D*,TH3D*,TH3D*,TH3D*> > p_hists;

    // loop over flavors
    for (Int_t f = elec; f <= nc; f++) {
      
      TH3D *h_NH_det_nu, *h_NH_det_nub, *h_IH_det_nu, *h_IH_det_nub;

      // set the pointers to the correct histograms
      if (f < nc) {
	
	TString hname_nu  = "detflux/detflux_" + fFlav_str[f] + "_cc_nu";
	TString hname_nub = "detflux/detflux_" + fFlav_str[f] + "_cc_nub";

	h_NH_det_nu  = (TH3D*)nhfile.Get(hname_nu);
	h_NH_det_nub = (TH3D*)nhfile.Get(hname_nub);
	h_IH_det_nu  = (TH3D*)ihfile.Get(hname_nu);
	h_IH_det_nub = (TH3D*)ihfile.Get(hname_nub);

      }

      else {

	h_NH_det_nu = (TH3D*)nhfile.Get("detflux/detflux_elec_nc_nu");
	h_NH_det_nu->Add( (TH3D*)nhfile.Get("detflux/detflux_muon_nc_nu") );
	h_NH_det_nu->Add( (TH3D*)nhfile.Get("detflux/detflux_tau_nc_nu") );

	h_NH_det_nub = (TH3D*)nhfile.Get("detflux/detflux_elec_nc_nub");
	h_NH_det_nub->Add( (TH3D*)nhfile.Get("detflux/detflux_muon_nc_nub") );
	h_NH_det_nub->Add( (TH3D*)nhfile.Get("detflux/detflux_tau_nc_nub") );

	h_IH_det_nu = (TH3D*)ihfile.Get("detflux/detflux_elec_nc_nu");
	h_IH_det_nu->Add( (TH3D*)ihfile.Get("detflux/detflux_muon_nc_nu") );
	h_IH_det_nu->Add( (TH3D*)ihfile.Get("detflux/detflux_tau_nc_nu") );

	h_IH_det_nub = (TH3D*)ihfile.Get("detflux/detflux_elec_nc_nub");
	h_IH_det_nub->Add( (TH3D*)ihfile.Get("detflux/detflux_muon_nc_nub") );
	h_IH_det_nub->Add( (TH3D*)ihfile.Get("detflux/detflux_tau_nc_nub") );

      }

      // per each flavor there are 4 histograms
      auto f_hists = std::make_tuple( (TH3D*)h_NH_det_nu->Clone(), (TH3D*)h_NH_det_nub->Clone(),
				      (TH3D*)h_IH_det_nu->Clone(), (TH3D*)h_IH_det_nub->Clone() );

      std::get<NH_nu> (f_hists)->SetDirectory(0);
      std::get<NH_nub>(f_hists)->SetDirectory(0);
      std::get<IH_nu> (f_hists)->SetDirectory(0);
      std::get<IH_nub>(f_hists)->SetDirectory(0);

      // add the histograms to the flavors vector
      p_hists.push_back( f_hists );
      
    } //end loop over flavors

    // add the flavors vector to the vector for pairs
    fFluxHists.push_back(p_hists);

  }

  cout << "NOTICE ASYM::InitFluxHists() finished initalising flux histograms." << endl;

}

//*****************************************************************************************

/**
   This function initialises the structure `ASYM::fSimHists` that are used for event weight calculation in filling the selections.
 */
void ASYM::InitVars() {

  // initialize weight histogram for each pair and flavor

  for (Int_t pair_index = 0; pair_index < fFluxHists.size(); pair_index++) {

    TH3D *h_template = std::get<NH_nu>(fFluxHists[pair_index][0]);
    vector< std::pair<TH3D*, TH3D*> > flav_hists;

    for (Int_t f = elec; f <= nc; f++) {

      TH3D* h_nu  = (TH3D*)h_template->Clone();
      TH3D* h_nub = (TH3D*)h_template->Clone();
      h_nu->Reset();
      h_nub->Reset();

      TString hname_nu  = "h_w_" + fFlav_str[f] + "_nu";
      TString hname_nub = "h_w_" + fFlav_str[f] + "_nub";
      h_nu ->SetNameTitle(hname_nu , hname_nu);
      h_nub->SetNameTitle(hname_nub, hname_nub);
      h_nu->SetDirectory(0);
      h_nub->SetDirectory(0);

      flav_hists.push_back( std::make_pair(h_nu, h_nub) );

    }

    fSimHists.push_back( flav_hists );

  }

  cout << "NOTICE ASYM::InitVars() finished initalising variables." << endl; 
 
}

//*****************************************************************************************

/**
   This function fills the `ASYM::fSimHists` for event weighting and the detector responses in `ASYM::fResponse`.
 */
void ASYM::FillSimHists() {

  for (Int_t i = 0; i < fSp->GetTree()->GetEntries(); i++) {
    
    fSp->GetTree()->GetEntry(i);
    SummaryEvent *evt = fSp->GetEvt();
    
    // ignore muons and noise
    if ( evt->Get_MC_is_neutrino() < 0.5 ) continue;

    // fill the detector responses
    for (auto &rvec: fResponse) {
      for (auto &r: rvec) r->Fill(evt);
    }

    // charged current
    if ( evt->Get_MC_is_CC() > 0.5 ) {

      // check for the true flavor of the neutrino
      for (Int_t f = elec; f <= tau; f++) {

	if ( TMath::Abs( evt->Get_MC_type() ) == fFlav_pdg[f] ) {

	  // loop simulation histograms for each NH_IH pair
	  for (auto &simh_pair: fSimHists) {

	    // fill hist, depending on whether it is nu/nub
	    if ( evt->Get_MC_type() > 0 ) {
	      simh_pair[f].first ->Fill( evt->Get_MC_energy(), -evt->Get_MC_dir_z(), evt->Get_MC_bjorkeny() );
	    }
	    else {
	      simh_pair[f].second->Fill( evt->Get_MC_energy(), -evt->Get_MC_dir_z(), evt->Get_MC_bjorkeny() );
	    }

	  }

	}

      } // end loop over flavors

    }
    // neutral current
    else {

      // loop simulation histograms for each NH_IH pair
      for (auto &simh_pair: fSimHists) {

	// fill nc hist, depending on whether it is nu/nub
	if ( evt->Get_MC_type() > 0 ) {
	  simh_pair[nc].first ->Fill( evt->Get_MC_energy(), -evt->Get_MC_dir_z(), evt->Get_MC_bjorkeny() );
	}
	else {
	  simh_pair[nc].second->Fill( evt->Get_MC_energy(), -evt->Get_MC_dir_z(), evt->Get_MC_bjorkeny() );
	}

      }

    }

  }

  cout << "NOTICE ASYM::FillSimHists() finished filling histograms for weighting." << endl;

}

//*****************************************************************************************

/**
   This function initialises the structure `ASYM::fEvtSels` and `ASYM::fResponse`.
 */
void ASYM::InitEvtSels_1() {

  for (Int_t i = 0; i < fFluxHists.size(); i++) {

    TH3D* h_template = std::get<NH_nu>(fFluxHists[i][0]);

    Double_t ebins  = h_template->GetXaxis()->GetNbins();
    Double_t emin   = h_template->GetXaxis()->GetBinLowEdge(1);
    Double_t emax   = h_template->GetXaxis()->GetBinLowEdge( ebins + 1 );
    
    Double_t ctbins = h_template->GetYaxis()->GetNbins();
    Double_t ctmin  = h_template->GetYaxis()->GetBinLowEdge(1);
    Double_t ctmax  = h_template->GetYaxis()->GetBinLowEdge( ctbins + 1 );
    
    Double_t bybins = h_template->GetZaxis()->GetNbins();
    Double_t bymin  = h_template->GetZaxis()->GetBinLowEdge(1);
    Double_t bymax  = h_template->GetZaxis()->GetBinLowEdge( bybins + 1 );

    // event selection template for tracks
    EventSelection t(EventSelection::track, "track", NULL, ebins, emin, emax, ctbins, ctmin, ctmax, bybins, bymin, bymax);
    t.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
    t.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
    t.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
    t.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
    t.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

    // detector response for tracks
    DetResponse r_t(DetResponse::track, "track", ebins, emin, emax, ctbins, ctmin, ctmax, bybins, bymin, bymax);
    r_t.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,  0.5, true );
    r_t.AddCut( &SummaryEvent::Get_track_ql1       , std::greater<double>()   ,  0.5, true );
    r_t.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   ,  0.6, true );
    r_t.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(), 0.05, true );
    r_t.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), 0.18, true );

    // event selection template for showers
    EventSelection s(EventSelection::shower, "shower", NULL, ebins, emin, emax, ctbins, ctmin, ctmax, bybins, bymin, bymax);
    s.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
    s.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
    s.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(),  0.6, true );
    s.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
    s.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

    // detector response for showers
    DetResponse r_s(DetResponse::shower, "shower", ebins, emin, emax, ctbins, ctmin, ctmax, bybins, bymin, bymax);
    r_s.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,  0.5, true );
    r_s.AddCut( &SummaryEvent::Get_shower_ql1     , std::greater<double>()   ,  0.5, true );
    r_s.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(),  0.6, true );
    r_s.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(), 0.05, true );
    r_s.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(),  0.5, true );

    vector<EventSelection*> NH_sels = { new EventSelection(t), new EventSelection(s) };
    vector<EventSelection*> IH_sels = { new EventSelection(t), new EventSelection(s) };
    vector<DetResponse*>    Resps   = { new DetResponse(r_t) , new DetResponse(r_s)  };
    fEvtSels.push_back( std::make_pair(NH_sels, IH_sels) );
    fResponse.push_back( Resps );

  }

  cout << "NOTICE ASYM::InitEvtSels_1() finished initalising event selections." << endl;

}

//*****************************************************************************************

/**
   This function fills the event selections `ASYM::fEvtSels`.
 */
void ASYM::FillSelections() {

  for (Int_t i = 0; i < fSp->GetTree()->GetEntries(); i++) {

    fSp->GetTree()->GetEntry(i);
    SummaryEvent *evt = fSp->GetEvt();

    if ( evt->Get_MC_is_neutrino() < 0.5 ) continue;
      
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
      TH3D *h_NH_det_nu  = std::get<NH_nu> (fFluxHists[p][flav]);
      TH3D *h_NH_det_nub = std::get<NH_nub>(fFluxHists[p][flav]);
      TH3D *h_IH_det_nu  = std::get<IH_nu> (fFluxHists[p][flav]);
      TH3D *h_IH_det_nub = std::get<IH_nub>(fFluxHists[p][flav]);

      TH3D* h_sim_nu  = fSimHists[p][flav].first;
      TH3D* h_sim_nub = fSimHists[p][flav].second;

      // get bins for weight calculation
      Int_t ebin  = h_NH_det_nu->GetXaxis()->FindBin(  evt->Get_MC_energy()   );
      Int_t ctbin = h_NH_det_nu->GetYaxis()->FindBin( -evt->Get_MC_dir_z()    );
      Int_t bybin = h_NH_det_nu->GetZaxis()->FindBin(  evt->Get_MC_bjorkeny() );

      // calculate weights for this event, NH and IH
      Double_t w_NH, w_IH;
      if ( evt->Get_MC_type() > 0 ) {
	w_NH = h_NH_det_nu->GetBinContent(ebin, ctbin, bybin)/h_sim_nu->GetBinContent(ebin, ctbin, bybin);
	w_IH = h_IH_det_nu->GetBinContent(ebin, ctbin, bybin)/h_sim_nu->GetBinContent(ebin, ctbin, bybin);
      }
      else {
	w_NH = h_NH_det_nub->GetBinContent(ebin, ctbin, bybin)/h_sim_nub->GetBinContent(ebin, ctbin, bybin);
	w_IH = h_IH_det_nub->GetBinContent(ebin, ctbin, bybin)/h_sim_nub->GetBinContent(ebin, ctbin, bybin);
      }

      // call fill for all event classes for NH and IH with this event
      for (auto ev: fEvtSels[p].first)  ev->Fill(evt, w_NH);
      for (auto ev: fEvtSels[p].second) ev->Fill(evt, w_IH);

    }

  }

  cout << "NOTICE ASYM::FillSelections() finished filling event selections." << endl;

}

//*****************************************************************************************

/**
   This function creates and fills the histograms `ASYM::fDetHists` using the detector responses stored in `ASYM::fResponse`.

   The histograms serve to cross-check that two methods - one using event weights calculated from `ASYM::fSimHists` and the other using the `DetResponse` class - give identical results.
 */
void ASYM::FillFromResponse() {

  for (Int_t pair_index = 0; pair_index < fEvtSels.size(); pair_index++) {

    if ( fEvtSels[pair_index].first.size() != fEvtSels[pair_index].second.size() ) {
      cout << "ERROR! ASYM::FillFromResponse() assumes the same number of event selections for NH/IH!" << endl;
      return;
    }
    
    Int_t Nevtsels = fEvtSels[pair_index].first.size();
    vector< TH3D* > hists_NH, hists_IH;

    // loop over the event selections
    for (Int_t evsel_index = 0; evsel_index < Nevtsels; evsel_index++) {

      // get the response
      DetResponse *detresp = fResponse[pair_index][evsel_index];
      
      // create a 3D histogram as in the event selection for NH and IH
      TH3D *evselh_NH  = fEvtSels[pair_index].first[evsel_index]->Get_h_E_costh_by();
      TH3D *evselh_IH  = fEvtSels[pair_index].second[evsel_index]->Get_h_E_costh_by();
      TString hname_NH = "NH3D_resp_" + fEvtSels[pair_index].first[evsel_index]->Get_SelName();
      TString hname_IH = "IH3D_resp_" + fEvtSels[pair_index].second[evsel_index]->Get_SelName();
      TH3D *resph_NH   = (TH3D*)evselh_NH->Clone();
      TH3D *resph_IH   = (TH3D*)evselh_IH->Clone();
      resph_NH->SetNameTitle(hname_NH, hname_NH);
      resph_IH->SetNameTitle(hname_IH, hname_IH);
      resph_NH->Reset();
      resph_IH->Reset();

      // fill the 'detected' histograms for NH and IH using the response
      for (Int_t ebin = 1; ebin <= resph_NH->GetXaxis()->GetNbins(); ebin++) {
	for (Int_t ctbin = 1; ctbin <= resph_NH->GetYaxis()->GetNbins(); ctbin++) {
	  for (Int_t bybin = 1; bybin <= resph_NH->GetZaxis()->GetNbins(); bybin++) {
	    
	    Double_t E  = resph_NH->GetXaxis()->GetBinCenter(ebin);
	    Double_t ct = resph_NH->GetYaxis()->GetBinCenter(ctbin);
	    Double_t by = resph_NH->GetZaxis()->GetBinCenter(bybin);

	    auto true_bins = detresp->GetBinWeights(E, ct, by);

	    Double_t bc_NH = 0;
	    Double_t bc_IH = 0;
	    Double_t bc_NH_err = 0;
	    Double_t bc_IH_err = 0;
	    for (auto &tb: true_bins) {
	      TH3D *detected_NH, *detected_IH;
	      Int_t flav = tb.fFlav;
	      if (tb.fIsCC == 0) flav = nc;

	      if (tb.fIsNB == 0) { 
		detected_NH = std::get<NH_nu> (fFluxHists[pair_index][flav]); 
		detected_IH = std::get<IH_nu> (fFluxHists[pair_index][flav]); 
	      }
	      else { 
		detected_NH = std::get<NH_nub> (fFluxHists[pair_index][flav]); 
		detected_IH = std::get<IH_nub> (fFluxHists[pair_index][flav]); 
	      }

	      bc_NH += tb.fW * detected_NH->GetBinContent(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin);
	      bc_IH += tb.fW * detected_IH->GetBinContent(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin);

	      bc_NH_err += TMath::Power( tb.fWE * detected_NH->GetBinContent(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin), 2 );
	      bc_IH_err += TMath::Power( tb.fWE * detected_IH->GetBinContent(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin), 2 );
	    }

	    resph_NH->SetBinContent(ebin, ctbin, bybin, bc_NH);
	    resph_NH->SetBinError(ebin, ctbin, bybin, TMath::Sqrt(bc_NH_err));
	    resph_IH->SetBinContent(ebin, ctbin, bybin, bc_IH);
	    resph_IH->SetBinError(ebin, ctbin, bybin, TMath::Sqrt(bc_IH_err));

	  } // end loop over z bins
	} // end loop over y bins
      } // end loop over x bins

      hists_NH.push_back(resph_NH);
      hists_IH.push_back(resph_IH);

    } // end loop over event selections

    fDetHists.push_back( std::make_pair(hists_NH, hists_IH) );

  } // end loop over fEvtSels. fEvtSels.first is a vector of event selections for NH file, .second for IH file

}

//*****************************************************************************************

/**
   This function calculates the asymmetries between `EventSelection`'s for NH and IH.
   \param NH_IH_pairs     Vector of filename pairs, vec[index].first is the NH file, vec[index].second the IH file.
   \param idstr           Indentifier added to the outputs.
 */
void ASYM::GetAsymmetries(vector< std::pair<TString, TString> > NH_IH_pairs, TString idstr) {

  // loop over NH/IH pairs
  for (UInt_t p = 0; p < fEvtSels.size(); p++) {

    vector<EventSelection*> NHsels = fEvtSels[p].first;
    vector<EventSelection*> IHsels = fEvtSels[p].second;
    vector<TH3D*> NHhists = fDetHists[p].first;
    vector<TH3D*> IHhists = fDetHists[p].second;

    if ( ( NHsels.size() != IHsels.size() ) || ( NHhists.size() != IHhists.size() ) || 
	 ( NHsels.size() != NHhists.size() ) ) {
      throw std::invalid_argument( " ERROR! ASYM::GetAsymmetries() unexpected event selections." );
    }

    // store the input FluxChain.C output filenames in the header
    FileHeader Asym("Asymmetry");
    Asym.AddParameter("NHflux", NH_IH_pairs[p].first);
    Asym.AddParameter("IHflux", NH_IH_pairs[p].second);
    vector<TObject*> hout;

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
      hout.push_back( NHsels[c]->Get_h_E_costh()->Clone("NH_" + selname) );
      hout.push_back( NHsels[c]->Get_h_E_costh_by()->Clone("NH3D_" + selname) );
      hout.push_back( NHhists[c]->Clone() );
      hout.push_back( IHsels[c]->Get_h_E_costh()->Clone("IH_" + selname) );
      hout.push_back( IHsels[c]->Get_h_E_costh_by()->Clone("IH3D_" + selname) );
      hout.push_back( IHhists[c]->Clone() );

    }

    // write the output to file
    TFile fout("output/Asymmetry_sel" + (TString)to_string(p) + "_" + idstr + ".root", "RECREATE");
    Asym.WriteHeader(&fout);
    for (auto &h: hout) { h->Write(); delete h; }
    fout.Close();

  }

  // write out the responses for investigation
  for (auto &rvec: fResponse) {
    for (auto &r: rvec) {
      TString name = "output/Response_by=" + (TString)to_string(r->GetHist3D()->GetZaxis()->GetNbins()) + r->Get_RespName() + ".root";
      r->WriteToFile(name);
    }
  }

  cout << "NOTICE ASYM::GetAsymmetries() finished calculating asymmetries." << endl;

}

//*****************************************************************************************

void ASYM::CleanUp() {

  if (fSp) delete fSp;
  
  for (auto &p: fSimHists) {
    for (Int_t f = 0; f < nc; f++) {
      if (p[f].first) delete p[f].first;
      if (p[f].second) delete p[f].second;
    }
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

  for (auto &rvec: fResponse) {
    for (auto &r: rvec) if (r) delete r;
  }

  cout << "NOTICE ASYM::CleanUp() finished." << endl;
}
