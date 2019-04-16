#include "DetResponse.h"
#include "NMHUtils.h"
#include "FileHeader.h"

#include "TAxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"

#include<iostream>
#include <stdexcept>
using namespace std;

/** Constructor.

    See `AbsResponse` constructor with matching parameter interface for more info.
    
 */
DetResponse::DetResponse(reco reco_type, TString resp_name,
			 Int_t ebins , Double_t emin , Double_t emax  ,
			 Int_t ctbins, Double_t ctmin, Double_t ctmax ,
			 Int_t bybins, Double_t bymin, Double_t bymax ) : AbsResponse(reco_type, resp_name, ebins, emin, emax, ctbins, ctmin, ctmax, bybins, bymin, bymax) {
  
  fNormalised = false;
    
  //----------------------------------------------------------
  // initialize histograms
  //----------------------------------------------------------
  
  for (auto &f: fFlavs) {
    for (auto &i: fInts) {
      for (auto &p: fPols) {
	TString hname = "hsim_" + f.second + "_" + i.second + "_" + p.second + "_" + fRespName;
	fhSim[f.first][i.first][p.first] = CloneFromTemplate( fhBinsTrue, hname );
      }
    }
  }

  fHResp         = CloneFromTemplate(fhBinsReco, "hresp_" + fRespName);
  fhAtmMuCount1y = CloneFromTemplate(fhBinsReco, "hAtmMuCount1y_" + fRespName);
  fhNoiseCount1y = CloneFromTemplate(fhBinsReco, "hNoiseCount1y_" + fRespName);
  
  //----------------------------------------------------------
  // calculate the binning for the fResp structure and init fResp
  // bin [0] is underflow, bin[ axis->GetNbins() ] is the last counting bin (hence dimension length +1),
  // bin[ axis->GetNbins()+1 ] is overflow (hence dimension length + 2)
  //----------------------------------------------------------
  fEbins  = fHResp->GetXaxis()->GetNbins() + 2;
  fCtbins = fHResp->GetYaxis()->GetNbins() + 2;
  fBybins = fHResp->GetZaxis()->GetNbins() + 2;

  InitResponse(fEbins, fCtbins, fBybins);
  
}

//*********************************************************************************

/**
   Copy constructor.
   \param name     name of the copy, should not match the name of the other, this may confuse ROOT histogram naming scheme
   \param detresp  the copied response instance
 */
DetResponse::DetResponse(TString name, const DetResponse &detresp) : AbsResponse(name, detresp) {

  fNormalised = detresp.fNormalised;
  fEbins      = detresp.fEbins;
  fCtbins     = detresp.fCtbins;
  fBybins     = detresp.fBybins;

  for (auto f: fFlavs) {
    for (auto i: fInts) {
      for (auto p: fPols) {
	TString hname = "hsim_" + f.second + "_" + i.second + "_" + p.second + "_" + fRespName;
	fhSim[f.first][i.first][p.first] = (TH3D*)detresp.fhSim[f.first][i.first][p.first]->Clone(hname);
      }
    }
  }

  fHResp = (TH3D*)detresp.fHResp->Clone("hresp_" + fRespName);
  fhAtmMuCount1y = (TH3D*)detresp.fhAtmMuCount1y->Clone("hAtmMuCount1y_" + fRespName);
  fhNoiseCount1y = (TH3D*)detresp.fhNoiseCount1y->Clone("hNoiseCount1y_" + fRespName);

  InitResponse(fEbins, fCtbins, fBybins);

  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {
      for (Int_t bybin = 0; bybin < fBybins; bybin++) {
	fResp[ebin][ctbin][bybin] = detresp.fResp[ebin][ctbin][bybin];
      }
    }
  }

}

//*********************************************************************************

/** Destructor. */
DetResponse::~DetResponse() {
  
  CleanResponse();

  if (fHResp) delete fHResp;
  if (fhAtmMuCount1y) delete fhAtmMuCount1y;
  if (fhNoiseCount1y) delete fhNoiseCount1y;

  for (auto &f: fFlavs) {
    for (auto &i: fInts) {
      for (auto &p: fPols) {
  	delete fhSim[f.first][i.first][p.first];
      }
    }
  }
  
  TIter next(&fHeap);
  TObject *obj = NULL;
  while ( (obj = next() ) ) if (obj) delete obj;

}

//*********************************************************************************

/** 
    Private function to create a histogram on the stack from template
    \param tmpl  Template histogram
    \param name  Name of the created histogram
    \return      Pointer to the new histogram
*/
TH3D* DetResponse::CloneFromTemplate(TH3D* tmpl, TString name) {

  TH3D* ret = (TH3D*)tmpl->Clone(name);
  ret->SetNameTitle(name, name);
  ret->Reset();
  ret->SetDirectory(0);

  return ret;

}

//*********************************************************************************

/**
   Private function to initialise the member `fResp` structure.
   \param ebins  Number of energy bins
   \param ctbins Number of cos-theta bins
   \param bybins Number of bjorken-y bins
 */
void DetResponse::InitResponse(Int_t ebins, Int_t ctbins, Int_t bybins) {

  fResp = new vector<TrueB>** [ebins]();
  for (Int_t ebin = 0; ebin < ebins; ebin++) {
    fResp[ebin] = new vector<TrueB>* [ctbins]();
    for (Int_t ctbin = 0; ctbin < ctbins; ctbin++) {
      fResp[ebin][ctbin] = new vector<TrueB> [bybins]();
    }
  }

}

//*********************************************************************************

/**
   Private function to de-allocate memory of the member `fResp` structure
 */
void DetResponse::CleanResponse() {

  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {
      if (fResp[ebin][ctbin]) delete[] fResp[ebin][ctbin];
    }
    if (fResp[ebin]) delete[] fResp[ebin];
  }
  delete[] fResp;

}

//*********************************************************************************

/**
   This function fills a `SummaryEvent` to the response structure.
   \param evt    Pointer to a `SummaryEvent`
 */
void DetResponse::Fill(SummaryEvent *evt) {

  if (fNormalised) {
    throw std::invalid_argument( "ERROR! DetResponse::Fill() cannot fill an already normalised response!" );
  }
  
  // determine event type (either neutrinos or other). Catch the exception if the particle type is
  // not in the `fType_to_Supported` map and throw `invalid_argument`.
  UInt_t flav;
  try {
    flav  = fType_to_Supported.at( (UInt_t)TMath::Abs( evt->Get_MC_type() ) );
  }
  catch (const std::out_of_range& oor) {
    throw std::invalid_argument("ERROR! DetResponse::Fill() unknown particle type " +
				to_string( TMath::Abs( evt->Get_MC_type() ) ) );
  }
  UInt_t is_cc = evt->Get_MC_is_CC();
  UInt_t is_nb = (UInt_t)(evt->Get_MC_type() < 0 );


  // logic for dealing with neutrino events
  if (flav <= TAU) {
    FillNuEvents(flav, is_cc, is_nb, evt);
  }
  // logic for dealing with other (atmospheric muons, noise) events
  else {

    // if event does not pass the cuts return
    if ( !PassesCuts(evt) ) return;

    // set the reconstruction observables
    SetObservables(evt); //implemented in EventFilter.C
  
    if      ( flav == ATMMU ) fhAtmMuCount1y->Fill( fEnergy, -fDir.z(), fBy, evt->Get_MC_w1y() );
    else if ( flav == NOISE ) fhNoiseCount1y->Fill( fEnergy, -fDir.z(), fBy, evt->Get_MC_w1y() );
    else {
      throw std::invalid_argument( "ERROR! DetResponse::Fill() unknown particle with flavor " + to_string(flav) );
    }

  }
 
}

//*********************************************************************************

/**
   Internal function to fill neutrino `SummaryEvent`'s to the response structure.
   \param flav   Neutrino flavor
   \param is_cc  cc interaction flag
   \param is_nb  is nu-bar flag
   \param evt    Pointer to summary event
 */
void DetResponse::FillNuEvents(UInt_t flav, UInt_t is_cc, UInt_t is_nb, SummaryEvent *evt) {

  //-----------------------------------------------------------------------
  // exclude events that are outside of the simulation range
  //-----------------------------------------------------------------------
  Double_t nebins  = fhSim[flav][is_cc][is_nb]->GetXaxis()->GetNbins();
  Double_t emin    = fhSim[flav][is_cc][is_nb]->GetXaxis()->GetBinLowEdge(1);
  Double_t emax    = fhSim[flav][is_cc][is_nb]->GetXaxis()->GetBinUpEdge(nebins);

  Double_t nctbins = fhSim[flav][is_cc][is_nb]->GetYaxis()->GetNbins();
  Double_t ctmin   = fhSim[flav][is_cc][is_nb]->GetYaxis()->GetBinLowEdge(1);
  Double_t ctmax   = fhSim[flav][is_cc][is_nb]->GetYaxis()->GetBinUpEdge(nctbins);

  Double_t nbybins = fhSim[flav][is_cc][is_nb]->GetZaxis()->GetNbins();
  Double_t bymin   = fhSim[flav][is_cc][is_nb]->GetZaxis()->GetBinLowEdge(1);
  Double_t bymax   = fhSim[flav][is_cc][is_nb]->GetZaxis()->GetBinUpEdge(nbybins);

  if (  evt->Get_MC_energy()   < emin  ||  evt->Get_MC_energy()   >= emax  ||
       -evt->Get_MC_dir_z()    < ctmin || -evt->Get_MC_dir_z()    >= ctmax ||
	evt->Get_MC_bjorkeny() < bymin ||  evt->Get_MC_bjorkeny() >= bymax ) {
    return;
  }

  //-----------------------------------------------------------------------
  // fill the histogram that counts the total number of simulated events
  //-----------------------------------------------------------------------
  fhSim[flav][is_cc][is_nb]->Fill( evt->Get_MC_energy(), -evt->Get_MC_dir_z(), evt->Get_MC_bjorkeny() );

  //-----------------------------------------------------------------------
  // if the event does not pass cuts return; otherwise set the observables depending
  // on the selected reco type of this class and fill the detector response
  //-----------------------------------------------------------------------
  if ( !PassesCuts(evt) ) return;

  SetObservables(evt); //implemented in EventFilter.C

  fHResp->Fill(fEnergy, -fDir.z(), fBy);

  Int_t  e_true_bin = fHResp->GetXaxis()->FindBin(  evt->Get_MC_energy()   );
  Int_t ct_true_bin = fHResp->GetYaxis()->FindBin( -evt->Get_MC_dir_z()    );
  Int_t by_true_bin = fHResp->GetZaxis()->FindBin(  evt->Get_MC_bjorkeny() );

  Int_t  e_reco_bin = fHResp->GetXaxis()->FindBin(  fEnergy   );
  Int_t ct_reco_bin = fHResp->GetYaxis()->FindBin( -fDir.z()  );
  Int_t by_reco_bin = fHResp->GetZaxis()->FindBin(  fBy       );

  TrueB tbin(flav, is_cc, is_nb, e_true_bin, ct_true_bin, by_true_bin);
  std::vector<TrueB>& true_bins = fResp[e_reco_bin][ct_reco_bin][by_reco_bin];
  auto index = std::find( true_bins.begin(), true_bins.end(), tbin );

  if ( index == true_bins.end() ) true_bins.push_back(tbin);
  else index->Increment();

}

//*********************************************************************************

/**
   This function normalises the response and needs to be called after all `SummaryEvent`'s have been filled to the response.
 */
void DetResponse::Normalise() {

  if (fNormalised) {
    cout << "WARNING! DetResponse::Normalise() response already normalised, exiting function." << endl;
    return;
  }
  
  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {
      for (Int_t bybin = 0; bybin < fBybins; bybin++) {
	
	std::vector<TrueB>& true_bins = fResp[ebin][ctbin][bybin];

	for (auto &tb: true_bins) {

	  //------------------------------------------------------------
	  // calculate the weight (fW) and weight MC error (fWE)
	  //------------------------------------------------------------

	  // histogram with total number of simulated events in E_true, ct_true, by_true bins
	  TH3D *hsim = fhSim[tb.fFlav][tb.fIsCC][tb.fIsNB];

	  // total number of simulated events in E_true, ct_true, by_true bin
	  Double_t nsim  = hsim->GetBinContent(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin);

	  // number of events that entered from true bin to reco bin (stored in tb.fW at this point)
	  Double_t nsel = tb.fW;
	  
	  // calculate the weight (fraction of simulated events in true bin that entered the reco bin)
	  tb.fW = nsel/nsim;
	  
	  // calculate the weight error from error propagation for fW = nsel/nsim under the
	  // assumption that nsel << nsim
	  tb.fWE = TMath::Sqrt(nsel)/nsim;
	  
	}
	
      }
    }
  }

  fNormalised = true; //set the flag to indicate that the response is normalised
}

//*********************************************************************************

/**
   Returns a vector of `TrueB`'s (true e, cos-theta, by-bins) that contribute to this e_reco, ct_reco, by_reco bin for neutrino events. 
   
   \param E_reco    Reconstructed energy
   \param ct_reco   Reconstructed cos-theta
   \param by_reco   Reconstructed bjorken-y
   \return          vector of `TrueB`'s that contribute to the reco bin.
 */
const std::vector<TrueB>& DetResponse::GetBinWeights(Double_t E_reco, Double_t ct_reco, Double_t by_reco) {

  if (!fNormalised) Normalise();
  
  Int_t ebin  = fHResp->GetXaxis()->FindBin(E_reco);
  Int_t ctbin = fHResp->GetYaxis()->FindBin(ct_reco);
  Int_t bybin = fHResp->GetZaxis()->FindBin(by_reco);

  return fResp[ebin][ctbin][bybin];
}

//*********************************************************************************

/**
   Returns a vector of `TrueB`'s (true e, cos-theta, by-bins) that contribute to this event for neutrino events.
   
   \param evt       Pointer to a summary event
   \return          vector of `TrueB`'s that contribute to the reco bin of this event.
 */
const std::vector<TrueB>& DetResponse::GetBinWeights(SummaryEvent *evt) {

  SetObservables(evt);
  return GetBinWeights( fEnergy, -fDir.z(), fBy );
  
}

//*********************************************************************************

/**
   Returns the number of atmospheric muon events in 1 year in the bin specified by the reconstruction variables.
   
   \param E_reco    Reconstructed energy
   \param ct_reco   Reconstructed cos-theta
   \param by_reco   Reconstructed bjorken-y
   \return          a pair with the atmospheric muon count in 1 year (first) and the MC statistical error (second)
 */
std::pair<Double_t, Double_t> DetResponse::GetAtmMuCount1y(Double_t E_reco, Double_t ct_reco, Double_t by_reco) {

  Int_t ebin  = fhAtmMuCount1y->GetXaxis()->FindBin(E_reco);
  Int_t ctbin = fhAtmMuCount1y->GetYaxis()->FindBin(ct_reco);
  Int_t bybin = fhAtmMuCount1y->GetZaxis()->FindBin(by_reco);
  
  return std::make_pair( fhAtmMuCount1y->GetBinContent(ebin, ctbin, bybin), fhAtmMuCount1y->GetBinError(ebin, ctbin, bybin) );

}

//*********************************************************************************

/**
   Returns the number of noise events in 1 year in the bin specified by the reconstruction variables.
   
   \param E_reco    Reconstructed energy
   \param ct_reco   Reconstructed cos-theta
   \param by_reco   Reconstructed bjorken-y
   \return          a pair with noise count in 1 year (first) and the MC statistical error (second)
 */
std::pair<Double_t, Double_t> DetResponse::GetNoiseCount1y(Double_t E_reco, Double_t ct_reco, Double_t by_reco) {

  Int_t ebin  = fhNoiseCount1y->GetXaxis()->FindBin(E_reco);
  Int_t ctbin = fhNoiseCount1y->GetYaxis()->FindBin(ct_reco);
  Int_t bybin = fhNoiseCount1y->GetZaxis()->FindBin(by_reco);
  
  return std::make_pair( fhNoiseCount1y->GetBinContent(ebin, ctbin, bybin), fhNoiseCount1y->GetBinError(ebin, ctbin, bybin) );

}

//*********************************************************************************

/**
   Function to write the response to file.

   \param filename  Name of the output file.

 */
void DetResponse::WriteToFile(TString filename) {

  if (!fNormalised) Normalise();

  FileHeader h("DetResponse");
  h.AddParameter( "fRespName", fRespName);
  h.AddParameter( "fEbins"   , (TString)to_string(fEbins) );
  h.AddParameter( "fCtbins"  , (TString)to_string(fCtbins) );
  h.AddParameter( "fBybins"  , (TString)to_string(fBybins) );
  
  TFile fout(filename, "RECREATE");
  TTree tout("detresponse","Detector response data");

  TrueB *tb = new TrueB();
  Int_t E_reco_bin, ct_reco_bin, by_reco_bin;

  tout.Branch("E_reco_bin" , &E_reco_bin , "E_reco_bin/I");
  tout.Branch("ct_reco_bin", &ct_reco_bin, "ct_reco_bin/I");
  tout.Branch("by_reco_bin", &by_reco_bin, "by_reco_bin/I");
  tout.Branch("TrueB", &tb, 2);

  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {
      for (Int_t bybin = 0; bybin < fBybins; bybin++) {

	E_reco_bin  = ebin;
	ct_reco_bin = ctbin;
	by_reco_bin = bybin;

	for (auto true_bin: fResp[ebin][ctbin][bybin]) {
	  *tb = true_bin;
	  tout.Fill();
	}

      }
    }
  }

  tout.Write();

  for (auto f: fFlavs) {
    for (auto i: fInts) {
      for (auto p: fPols) {
	// write without fRespName for easier read-in
	TString hname = f.second + "_" + i.second + "_" + p.second;
	fhSim[f.first][i.first][p.first]->Write(hname);
      }
    }
  }

  fHResp->Write("hresp");
  fhAtmMuCount1y->Write("atmmucount");
  fhNoiseCount1y->Write("noisecount");

  h.WriteHeader(&fout);

  fout.Close();
  delete tb;
  
}

//*********************************************************************************

/**
   Function to read the response from file.

   \param filename   Name of the file where the response is stored.
 */
void DetResponse::ReadFromFile(TString filename) {

  // if the detresponse is read from file, the binning is changed to match that
  // of the detresponse that is being read in
  CleanResponse();

  FileHeader h("for_read_in");
  h.ReadHeader(filename);

  fRespName   = h.GetParameter("fRespName");
  fNormalised = true;
  fEbins      = std::stoi( (string)h.GetParameter("fEbins")  );
  fCtbins     = std::stoi( (string)h.GetParameter("fCtbins") );
  fBybins     = std::stoi( (string)h.GetParameter("fBybins") );

  InitResponse(fEbins, fCtbins, fBybins);

  // read in the true bin data to the response

  TFile fin(filename, "READ");
  TTree *tin = (TTree*)fin.Get("detresponse");

  TrueB *tb = new TrueB();
  Int_t E_reco_bin, ct_reco_bin, by_reco_bin;

  tin->SetBranchAddress("E_reco_bin", &E_reco_bin);
  tin->SetBranchAddress("ct_reco_bin", &ct_reco_bin);
  tin->SetBranchAddress("by_reco_bin", &by_reco_bin);
  tin->SetBranchAddress("TrueB", &tb);

  for (Int_t i = 0; i < tin->GetEntries(); i++) {
    tin->GetEntry(i);
    fResp[E_reco_bin][ct_reco_bin][by_reco_bin].push_back( TrueB(*tb) );
  }

  // get the histgrams

  for (auto f: fFlavs) {
    for (auto i: fInts) {
      for (auto p: fPols) {
	TString hname1 = f.second + "_" + i.second + "_" + p.second; //name in file
	TString hname2 = "hsim_" + f.second + "_" + i.second + "_" + p.second + "_" + fRespName; //name for instance
	fhSim[f.first][i.first][p.first] = (TH3D*)fin.Get(hname1)->Clone();
	fhSim[f.first][i.first][p.first]->SetDirectory(0);
	fhSim[f.first][i.first][p.first]->SetNameTitle(hname2, hname2);
      }
    }
  }

  fHResp = (TH3D*)fin.Get("hresp")->Clone();
  fHResp->SetDirectory(0);
  fHResp->SetNameTitle("hresp_" + fRespName, "hresp_" + fRespName);

  fhAtmMuCount1y = (TH3D*)fin.Get("atmmucount")->Clone();
  fhAtmMuCount1y->SetDirectory(0);
  fhAtmMuCount1y->SetNameTitle("hAtmMuCount1y_" + fRespName, "hAtmMuCount1y_" + fRespName);

  fhNoiseCount1y = (TH3D*)fin.Get("noisecount")->Clone();
  fhNoiseCount1y->SetDirectory(0);
  fhNoiseCount1y->SetNameTitle("hNoiseCount1y_" + fRespName, "hNoiseCount1y_" + fRespName);

  delete tb;
  fin.Close();

}

//*********************************************************************************

/**
   Function to visualise the response in 2D.

   This function shows the weights of the true bins that contribute to the reco bin in 2D and the number of events in true bins that contribute to the reco bin in 2D. As the response itself is in 3D, summation over bjorken-y is performed, which makes the structure of the function somewhat cumbersome.

   \param e_reco    Reconstructed energy
   \param ct_reco   Reconstructed cos-theta
   \param outname   If specified, the histograms and canvases are written to output
   \return          A tuple with pointers to three canvases where the visualization is drawn.
 */
std::tuple<TCanvas*,TCanvas*,TCanvas*> DetResponse::DisplayResponse(Double_t e_reco, Double_t ct_reco, TString outname) {

  if (!fNormalised) Normalise();

  //----------------------------------------------------------------------
  // initialise the histograms for displaying the response
  //----------------------------------------------------------------------
  TH2D *hw[fFlavs.size()][fInts.size()][fPols.size()]; //histogram displaying weights by nu type
  TH2D *hn[fFlavs.size()][fInts.size()][fPols.size()]; //histogram displaying event counts by nu type

  TString suffix = "_Ereco=" + (TString)to_string(e_reco) + "_ctreco=" + (TString)to_string(ct_reco);
  
  for (auto f: fFlavs) {
    for (auto i: fInts) {
      for (auto p: fPols) {
	TString wname = "weights_" + f.second + "_" + i.second + "_" + p.second + suffix;
	TString nname = "counts_" + f.second + "_" + i.second + "_" + p.second + suffix;
	hw[f.first][i.first][p.first] = (TH2D*)fHResp->Project3D("yx")->Clone();
	hn[f.first][i.first][p.first] = (TH2D*)fHResp->Project3D("yx")->Clone();
	hw[f.first][i.first][p.first]->Reset();
	hn[f.first][i.first][p.first]->Reset();
	hw[f.first][i.first][p.first]->SetDirectory(0);
	hn[f.first][i.first][p.first]->SetDirectory(0);
	hw[f.first][i.first][p.first]->SetNameTitle(wname,wname);
	hn[f.first][i.first][p.first]->SetNameTitle(nname,nname);
	fHeap.Add(hw[f.first][i.first][p.first]);
	fHeap.Add(hn[f.first][i.first][p.first]);
      }
    }
  }

  //----------------------------------------------------------------------
  // histogram displaying the number of events in the selected reconstruction bin
  //----------------------------------------------------------------------  
  TH2D *h_rb  = (TH2D*)fHResp->Project3D("yx")->Clone();
  h_rb->SetDirectory(0);
  h_rb->SetNameTitle("reco_bin"+suffix,"reco_bin"+suffix);
  Int_t ebin  = h_rb->GetXaxis()->FindBin(e_reco);
  Int_t ctbin = h_rb->GetYaxis()->FindBin(ct_reco);
  Double_t bc = h_rb->GetBinContent(ebin, ctbin);
  h_rb->Reset();
  h_rb->SetBinContent( ebin, ctbin, bc );

  //----------------------------------------------------------------------
  // sum over bjorken-y by adding all true bins over bjorken-y to the true_bins vector
  //----------------------------------------------------------------------
  vector<TrueB> true_bins;
  for (Int_t bybin = 0; bybin < fBybins; bybin++) {
    vector<TrueB>& tbs = fResp[ebin][ctbin][bybin];
    true_bins.insert( true_bins.end(), tbs.begin(), tbs.end() );
  }

  //----------------------------------------------------------------------
  // loop over true bins and add contents to the histograms visualizing the response
  //----------------------------------------------------------------------
  for (auto &tb: true_bins) {

    // total number of simulated events in the true bin is stored in fhSim
    Double_t totsim = fhSim[tb.fFlav][tb.fIsCC][tb.fIsNB]->GetBinContent( tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin );

    // number of events from true bin in reco bin - reconstructed from the weight
    Double_t recsim = totsim * tb.fW;

    // histogram that will display the weights in 2D
    hw[tb.fFlav][tb.fIsCC][tb.fIsNB]->SetBinContent( tb.fE_true_bin, tb.fCt_true_bin,
						     hw[tb.fFlav][tb.fIsCC][tb.fIsNB]->GetBinContent(tb.fE_true_bin, tb.fCt_true_bin) + recsim );

    // histogram that will display the number of event from true bin that entered the reco bin
    hn[tb.fFlav][tb.fIsCC][tb.fIsNB]->SetBinContent( tb.fE_true_bin, tb.fCt_true_bin,
						     hn[tb.fFlav][tb.fIsCC][tb.fIsNB]->GetBinContent(tb.fE_true_bin, tb.fCt_true_bin) + recsim );
  }

  //----------------------------------------------------------------------
  // now do the drawing; if outname is specified write histograms to output
  //----------------------------------------------------------------------

  TFile *fout = NULL;

  if (outname != "") {
    fout = new TFile(outname, "RECREATE");
  }
  
  Int_t padcount = 1;
  Int_t rows     = 3;
  Int_t cols     = 4;
  
  TCanvas *c1 = new TCanvas("Weights"+suffix,"Weights"+suffix, 1000, 600);
  c1->Divide(cols,rows);

  TCanvas *c2 = new TCanvas("Counts"+suffix,"Counts"+suffix, 1000, 600);
  c2->Divide(cols,rows);
  
  for (auto f: fFlavs) {
    for (auto i: fInts) {
      for (auto p: fPols) {

	if (padcount > rows*cols) {
	  throw std::invalid_argument( "ERROR! DetResponse::DisplayResponse() too few pads!" );
	}
	
	// the weights displaying histogram shows only numbers at this stage, these need to
	// be converted back to weights for a 2D response
	TH2D *hsim2D = (TH2D*)fhSim[f.first][i.first][p.first]->Project3D("yx")->Clone("clone");
	hw[f.first][i.first][p.first]->Divide(hsim2D);
	if (hsim2D) delete hsim2D;

	c1->cd(padcount);
	c1->GetPad(padcount)->SetLogx();
	hw[f.first][i.first][p.first]->Draw("colz");

	c2->cd(padcount);
	c2->GetPad(padcount)->SetLogx();
	hn[f.first][i.first][p.first]->Draw("text");

	if (fout) {
	  hw[f.first][i.first][p.first]->Write();
	  hn[f.first][i.first][p.first]->Write();
	}
	
	padcount++;
      }
    }
  }  
    
  TCanvas *c3 = new TCanvas("RecoBin"+suffix,"RecoBin"+suffix, 1000, 600);
  c3->SetLogx();
  h_rb->Draw("colz");
  
  fHeap.Add(h_rb);
  fHeap.Add(c1);
  fHeap.Add(c2);
  fHeap.Add(c3);

  if (fout) {
    fout->Close();
    delete fout;
  }
  
  return std::make_tuple(c1,c2,c3);
}

//*********************************************************************************

/** Function to get the histogram that counts the number of simulated events.

    \param flav   Neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param iscc   True for CC events, False for NC events
    \param isnb   True for nubar, False for nu
    \return       Pointer to a `TH3D` that counts the simulated events for the specified neutrino type
*/
TH3D* DetResponse::GetSimHist(UInt_t flav, Bool_t iscc, Bool_t isnb) const {

  if (flav > TAU) {
    throw std::invalid_argument("ERROR! DetResponse::GetSimHist() unknown flavor " + to_string(flav));
  }

  return fhSim[flav][(UInt_t)iscc][(UInt_t)isnb];
}
