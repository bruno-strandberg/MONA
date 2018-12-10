#include "NMHUtils.h"
#include<vector>
using namespace std;

void Normalised_BY_dists(TString byhists_olist, TString output = "out.root") {

  //-------------------------------------------------
  // merge the files
  //-------------------------------------------------

  vector<TString> files = NMHUtils::ReadLines(byhists_olist);

  TString syscmd = "hadd " + output;
  for (auto &fname: files) {
    syscmd += " " + fname;
  }

  Int_t sysret = 0;
  if ( NMHUtils::FileExists(output) ) sysret = system("rm " + output);
  sysret = system(syscmd);

  //-------------------------------------------------
  // read the merged histograms to memory
  //-------------------------------------------------

  TFile fin(output, "READ");
  vector< TH2D* > hist_list;

  for ( auto key: *fin.GetListOfKeys() ) {
    TH2D *h = (TH2D*)fin.Get( key->GetName() );
    hist_list.push_back( (TH2D*)h->Clone() );
    hist_list.back()->SetDirectory(0);
  }

  fin.Close();

  //-------------------------------------------------
  // normalise each histogram of the output at each energy bin
  //-------------------------------------------------
  
  for (auto h: hist_list) {

    for (Int_t ebin = 1; ebin <= h->GetXaxis()->GetNbins(); ebin++) {

      //project along Y axis to get bjorken-y dist at a certain energy
      TH1D* proj = h->ProjectionY("proj", ebin, ebin);
      if (proj->Integral() == 0.) continue;

      for (Int_t bybin = 1; bybin <= h->GetYaxis()->GetNbins(); bybin++) {
	Double_t bc = h->GetBinContent(ebin, bybin)/proj->Integral();
	Double_t be = h->GetBinError(ebin, bybin)/proj->Integral();
	h->SetBinContent(ebin, bybin, bc);
	h->SetBinError(ebin, bybin, be);
      } // end loop over y axis

    } // end loop over x axis

  } // end loop over histograms

  //-------------------------------------------------
  // overwrite the output with normalised histograms, add header for info
  //-------------------------------------------------
  FileHeader h("Normalised_BY_dists");
  h.AddParameter("gSeaGen_files","/in2p3/km3net/mc/atm_neutrino/KM3NeT_ORCA_115_23m_9m/v1.0/gSeaGen");

  TFile fout(output, "RECREATE");
  for (auto h: hist_list) {
    h->Write();
  }

  h.WriteHeader(&fout);

  fout.Close();

}
