#include "DetResponse.h"
#include "NMHUtils.h"

#include "TAxis.h"

#include<iostream>
using namespace std;

/** Constructor.
    
 */
DetResponse::DetResponse(reco reco_type, TString resp_name,
			 Int_t ebins , Double_t emin , Double_t emax  ,
			 Int_t ctbins, Double_t ctmin, Double_t ctmax ,
			 Int_t bybins, Double_t bymin, Double_t bymax ) : EventFilter(reco_type) {
  
  fRespName = resp_name;

  //----------------------------------------------------------
  // initialize axes for histograms
  //----------------------------------------------------------
  
  vector<Double_t> e_edges = NMHUtils::GetLogBins(ebins, emin, emax);  

  TAxis ct_axis( ctbins, ctmin, ctmax );
  TAxis by_axis( bybins, bymin, bymax );

  vector<Double_t> ct_edges, by_edges;
  
  for (Int_t ctbin = 1; ctbin <= ct_axis.GetNbins() + 1; ctbin++) {
    ct_edges.push_back( ct_axis.GetBinLowEdge(ctbin) );
  }

  for (Int_t bybin = 1; bybin <= by_axis.GetNbins() + 1; bybin++) {
    by_edges.push_back( by_axis.GetBinLowEdge(bybin) );
  }
  
  //----------------------------------------------------------
  // initialize histograms
  //----------------------------------------------------------
  
  for (auto &f: fFlavs) {
    for (auto &i: fInts) {
      for (auto &p: fPols) {
	TString hname = "hsim_" + f.second + "_" + i.second + "_" + p.second;
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
  
}

//*********************************************************************************

void DetResponse::Fill(SummaryEvent *evt) {

  //-----------------------------------------------------------------------
  // determine event type and fill hists that count the number
  // of simulated events per E_true, costh_true, by_true bin
  //-----------------------------------------------------------------------
  
  Int_t flav  = fType_to_Flav[ (Int_t)evt->Get_MC_type() ];
  Int_t is_cc = evt->Get_MC_is_CC();
  Int_t is_nb = (Int_t)(evt->Get_MC_type() < 0 );

  // for atm muons and noise use only one histogram
  if ( evt->Get_MC_is_neutrino() ) {
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

void DetResponse::Normalise() {

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
  
}

//*********************************************************************************

std::vector<TrueB>& DetResponse::GetBinWeights(Double_t E_reco, Double_t ct_reco, Double_t by_reco) {

  Int_t ebin  = fRespH->GetXaxis()->FindBin(E_reco);
  Int_t ctbin = fRespH->GetYaxis()->FindBin(ct_reco);
  Int_t bybin = fRespH->GetZaxis()->FindBin(by_reco);

  return fResp[ebin][ctbin][bybin];
}
