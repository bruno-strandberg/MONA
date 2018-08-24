#include "DetResponse.h"
#include "NMHUtils.h"

#include "TAxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"

#include<iostream>
#include <stdexcept>
using namespace std;

/** Constructor.

    \param reco_type         Type of reco variables used to fill the response, e.g. DetResponse::track; See `EventFilter`
    \param resp_name         Name of the response
    \param ebins             Number of energy bins
    \param emin              Energy minimum
    \param emax              Energy maximum
    \param ctbins            Number of cos-theta bins
    \param ctmin             cos-theta minimum
    \param ctmax             cos-theta maximum
    \param bybins            Number of bjorken-y bins
    \param bymin             bjorken-y minimum
    \param bymax             bjorken-y maximum
    
 */
DetResponse::DetResponse(reco reco_type, TString resp_name,
			 Int_t ebins , Double_t emin , Double_t emax  ,
			 Int_t ctbins, Double_t ctmin, Double_t ctmax ,
			 Int_t bybins, Double_t bymin, Double_t bymax ) : EventFilter(reco_type) {
  
  fRespName   = resp_name;
  fNormalised = false;
  
  //----------------------------------------------------------
  // initialize axes for histograms
  //----------------------------------------------------------
  
  vector<Double_t> e_edges  = NMHUtils::GetLogBins(ebins, emin, emax);  
  vector<Double_t> ct_edges = NMHUtils::GetBins(ctbins, ctmin, ctmax);
  vector<Double_t> by_edges = NMHUtils::GetBins(bybins, bymin, bymax);
  
  //----------------------------------------------------------
  // initialize histograms
  //----------------------------------------------------------
  
  for (auto &f: fFlavs) {
    for (auto &i: fInts) {
      for (auto &p: fPols) {
	TString hname = "hsim_" + f.second + "_" + i.second + "_" + p.second + "_" + fRespName;
	fhSim[f.first][i.first][p.first] = new TH3D(hname, hname,
						    ebins , &e_edges[0],
						    ctbins, &ct_edges[0],
						    bybins, &by_edges[0]);
      }
    }
  }

  fRespH = new TH3D("hresp_" + fRespName, "hresp_" + fRespName,
		    ebins , &e_edges[0],
		    ctbins, &ct_edges[0],
		    bybins, &by_edges[0]);

  //----------------------------------------------------------
  // calculate the binning for the fResp structure and init fResp
  // bin [0] is underflow, bin[ axis->GetNbins() ] is the last counting bin (hence dimension length +1),
  // bin[ axis->GetNbins()+1 ] is overflow (hence dimension length + 2)
  //----------------------------------------------------------
  fEbins  = fRespH->GetXaxis()->GetNbins() + 2;
  fCtbins = fRespH->GetYaxis()->GetNbins() + 2;
  fBybins = fRespH->GetZaxis()->GetNbins() + 2;

  fResp = new vector<TrueB>** [fEbins]();
  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    fResp[ebin] = new vector<TrueB>* [fCtbins]();
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {
      fResp[ebin][ctbin] = new vector<TrueB> [fBybins]();
    }
  }
  
}

//*********************************************************************************

/**
   Copy constructor.

   Without explicit cloning of the root histograms the destruction causes a seg fault.
 */
DetResponse::DetResponse(const DetResponse &detresp) : EventFilter(detresp) {

  fRespName   = detresp.fRespName;
  fNormalised = detresp.fNormalised;
  fEbins      = detresp.fEbins;
  fCtbins     = detresp.fCtbins;
  fBybins     = detresp.fBybins;

  for (auto f: fFlavs) {
    for (auto i: fInts) {
      for (auto p: fPols) {
	fhSim[f.first][i.first][p.first] = (TH3D*)detresp.fhSim[f.first][i.first][p.first]->Clone();
      }
    }
  }

  fRespH = (TH3D*)detresp.fRespH->Clone();
  
  fResp = new vector<TrueB>** [fEbins]();
  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    fResp[ebin] = new vector<TrueB>* [fCtbins]();
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {
      fResp[ebin][ctbin] = new vector<TrueB> [fBybins]();
      for (Int_t bybin = 0; bybin < fBybins; bybin++) {
	fResp[ebin][ctbin][bybin] = detresp.fResp[ebin][ctbin][bybin];
      }
    }
  }

}

//*********************************************************************************

/** Destructor. */
DetResponse::~DetResponse() {
  
  for (Int_t ebin = 0; ebin < fEbins; ebin++) {
    for (Int_t ctbin = 0; ctbin < fCtbins; ctbin++) {
      if (fResp[ebin][ctbin]) delete[] fResp[ebin][ctbin];
    }
    if (fResp[ebin]) delete[] fResp[ebin];
  }
  delete[] fResp;

  if (fRespH) delete fRespH;

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
   This function fills a `SummaryEvent` to the response structure.
   
   \param SummaryEvent    Pointer to a summary event
 */
void DetResponse::Fill(SummaryEvent *evt) {

  if (fNormalised) {
    throw std::invalid_argument( "ERROR! DetResponse::Fill() cannot fill an already normalised response!" );
  }
  
  //-----------------------------------------------------------------------
  // determine event type and fill hists that count the number
  // of simulated events per E_true, costh_true, by_true bin
  //-----------------------------------------------------------------------
  
  Int_t flav  = fType_to_Flav[ (Int_t)TMath::Abs( evt->Get_MC_type() ) ];
  Int_t is_cc = evt->Get_MC_is_CC();
  Int_t is_nb = (Int_t)(evt->Get_MC_type() < 0 );

  // for atm muons and noise use only one histogram
  if ( evt->Get_MC_is_neutrino() < 0.5 ) {
    is_cc = 0;
    is_nb = 0;
  }

  fhSim[flav][is_cc][is_nb]->Fill( evt->Get_MC_energy(), -evt->Get_MC_dir_z(), evt->Get_MC_bjorkeny() );

  //-----------------------------------------------------------------------
  // if the event does not pass cuts return; otherwise set the observables depending
  // on the selected reco type of this class and fill the detector response
  //-----------------------------------------------------------------------
  if ( !PassesCuts(evt) ) return;

  SetObservables(evt); //implemented in EventFilter.C

  fRespH->Fill(fEnergy, -fDir.z(), fBy);

  Int_t  e_true_bin = fRespH->GetXaxis()->FindBin(  evt->Get_MC_energy()   );
  Int_t ct_true_bin = fRespH->GetYaxis()->FindBin( -evt->Get_MC_dir_z()    );
  Int_t by_true_bin = fRespH->GetZaxis()->FindBin(  evt->Get_MC_bjorkeny() );

  Int_t  e_reco_bin = fRespH->GetXaxis()->FindBin(  fEnergy   );
  Int_t ct_reco_bin = fRespH->GetYaxis()->FindBin( -fDir.z()  );
  Int_t by_reco_bin = fRespH->GetZaxis()->FindBin(  fBy       );

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
	  TH3D *hsim = fhSim[tb.fFlav][tb.fIsCC][tb.fIsNB];
	  Double_t nsim  = hsim->GetBinContent(tb.fE_true_bin, tb.fCt_true_bin, tb.fBy_true_bin);
	  Double_t nreco = fRespH->GetBinContent(ebin, ctbin, bybin);
	  tb.fFracTrue = tb.fFracTrue/nsim;
	  tb.fFracReco = tb.fFracReco/nreco;
	}
	
      }
    }
  }

  fNormalised = true; //set the flag to indicate that the response is normalised
}

//*********************************************************************************

/**
   Returns a vector of `TrueB`'s (true e, cos-theta, by-bins) that contribute to this e_reco, ct_reco, by_reco bin. 

   NB! When parsing actual `SummaryEvent`'s, it is recommended to use `DetResponse::GetBinWeights(SummaryEvent *evt)`. This guarantees that the correct reco variables are used to look up the return vector.
   
   \param e_reco    Reconstructed energy
   \param t_reco    Reconstructed cos-theta
   \param by_reco   Reconstructed bjorken-y
   \return          vector of `TrueB`'s that contribute to the reco bin.
 */
std::vector<TrueB>& DetResponse::GetBinWeights(Double_t E_reco, Double_t ct_reco, Double_t by_reco) {

  if (!fNormalised) Normalise();
  
  Int_t ebin  = fRespH->GetXaxis()->FindBin(E_reco);
  Int_t ctbin = fRespH->GetYaxis()->FindBin(ct_reco);
  Int_t bybin = fRespH->GetZaxis()->FindBin(by_reco);

  return fResp[ebin][ctbin][bybin];
}

//*********************************************************************************

/**
   Returns a vector of `TrueB`'s (true e, cos-theta, by-bins) that contribute to this event.
   
   \param evt       Pointer to a summary event
   \return          vector of `TrueB`'s that contribute to the reco bin of this event.
 */
std::vector<TrueB>& DetResponse::GetBinWeights(SummaryEvent *evt) {

  SetObservables(evt);
  return GetBinWeights( fEnergy, -fDir.z(), fBy );
  
}

//*********************************************************************************

/**
   Function to write the response to file.

   \param filename  Name of the output file.

 */
void DetResponse::WriteToFile(TString filename) {

  if (!fNormalised) Normalise();
  
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

  fRespH->Write("hresp");

  fout.Close();
  delete tb;
  
}

//*********************************************************************************

/**
   Function to read the response from file.

   \param filename   Name of the file where the response is stored.
 */
void DetResponse::ReadFromFile(TString filename) {

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

  fRespH = (TH3D*)fin.Get("hresp")->Clone();
  fRespH->SetDirectory(0);
  fRespH->SetNameTitle("hresp_" + fRespName, "hresp_" + fRespName);

  delete tb;
  fin.Close();

  fNormalised = true;
}

//*********************************************************************************

/**
   Function to visualise the response in 2D.

   \param e_reco    Reconstructed energy
   \param ct_reco   Reconstructed bjorken-y
   \return          Pointer to a canvas where the visualisation is drawn.
 */
TCanvas* DetResponse::DisplayResponse(Double_t e_reco, Double_t ct_reco) {

  if (!fNormalised) Normalise();
  
  TH2D *h_rb  = (TH2D*)fRespH->Project3D("yx")->Clone();
  TH2D* h_tbs = (TH2D*)fRespH->Project3D("yx")->Clone();
  h_rb->Reset();
  h_tbs->Reset();
  h_rb->SetDirectory(0);
  h_tbs->SetDirectory(0);
  TString suffix = "_Ereco=" + (TString)to_string(e_reco) + "_ctreco=" + (TString)to_string(ct_reco);
  h_rb->SetNameTitle("reco_bin"+suffix,"reco_bin"+suffix);
  h_tbs->SetNameTitle("true_bins"+suffix,"true_bins"+suffix);

  Int_t ebin  = h_rb->GetXaxis()->FindBin(e_reco);
  Int_t ctbin = h_rb->GetYaxis()->FindBin(ct_reco);

  // sum over bjorken-y
  for (Int_t bybin = 0; bybin < fBybins; bybin++) {
        
    vector<TrueB>& true_bins = fResp[ebin][ctbin][bybin];
    if (true_bins.size() == 0) continue; //ignore underflow and overflow

    h_rb->SetBinContent( ebin, ctbin, h_rb->GetBinContent(ebin, ctbin) + 1 );

    for (auto &tb: true_bins) {
      h_tbs->SetBinContent( tb.fE_true_bin, tb.fCt_true_bin, 
			    h_tbs->GetBinContent(tb.fE_true_bin, tb.fCt_true_bin) + tb.fFracReco );
    }

  }

  TCanvas *c1 = new TCanvas("DispResponse"+suffix,"DispResponse"+suffix, 1000, 600);
  c1->Divide(2,1);
  c1->cd(1);
  c1->GetPad(1)->SetLogx();
  h_rb->Draw("colz");
  c1->cd(2);
  c1->GetPad(2)->SetLogx();
  h_tbs->Draw("colz");

  fHeap.Add(h_rb);
  fHeap.Add(h_tbs);
  fHeap.Add(c1);

  return c1;
}
