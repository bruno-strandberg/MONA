#include "NMHUtils.h"

/* This macro compared track reco energy resolution to shower reco energy resolution for muons*/
void resolutions(TString fname = "../../data/ORCA_MC_summary_ORCA7_23x9m_ECAP1018.root") {

  SummaryParser sp(fname);

  Int_t nebins = 10;
  
  auto ebins = NMHUtils::GetLogBins(nebins,1,100);
  
  TH2D *trkres = new TH2D("trkresolution","trkresolution", nebins, &ebins[0], 200,-100,100);
  TH2D *shwres = new TH2D("shwresolution","shwresolution", nebins, &ebins[0], 200,-100,100);

  TH2D *trkdir = new TH2D("trkdir","trkdir", 200, -1, 1, 200, -1, 1);
  TH2D *shwdir = new TH2D("shwdir","shwdir", 200, -1, 1, 200, -1, 1);

  TH2D *trkvstrue = new TH2D("trkvstrue","trkvstrue", 25,1,50,25,1,50);
  TH2D *shwvstrue = new TH2D("shwvstrue","shwvstrue", 25,1,50,25,1,50);

  
  // only look at muon-CC events where track reco "worked" (not NaN = track_ql0 > 0)
  EventFilter trackfilter;
  trackfilter.AddCut( &SummaryEvent::Get_MC_type   , std::equal_to<double>(),  14, false);
  trackfilter.AddCut( &SummaryEvent::Get_MC_type   , std::equal_to<double>(), -14, false);
  trackfilter.AddCut( &SummaryEvent::Get_track_ql0 , std::greater<double>() , 0.5, true);
  trackfilter.AddCut( &SummaryEvent::Get_MC_is_CC  , std::greater<double>() , 0.5, true);

  // only look at muon-CC events where shower reco "worked" (not NaN = shower_ql0 > 0)
  EventFilter showerfilter;
  showerfilter.AddCut( &SummaryEvent::Get_MC_type   , std::equal_to<double>(),  14, false);
  showerfilter.AddCut( &SummaryEvent::Get_MC_type   , std::equal_to<double>(), -14, false);
  showerfilter.AddCut( &SummaryEvent::Get_shower_ql0, std::greater<double>() , 0.5, true);
  showerfilter.AddCut( &SummaryEvent::Get_MC_is_CC  , std::greater<double>() , 0.5, true);

  // fill histograms
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    
    if (i%100000 == 0) cout << "Event: " << i << endl;

    SummaryEvent *evt = sp.GetEvt(i);

    if ( trackfilter.PassesCuts(evt) ) {
      trkdir->Fill( -evt->Get_MC_dir().z(), -evt->Get_track_dir().z() );
      trkvstrue->Fill( evt->Get_MC_energy(), evt->Get_track_energy() );
      trkres->Fill( evt->Get_MC_energy(), ( evt->Get_MC_energy() - evt->Get_track_energy() )/evt->Get_MC_energy() * 100 );
    }

    if ( showerfilter.PassesCuts(evt) ) {
      shwdir->Fill( -evt->Get_MC_dir().z(), -evt->Get_shower_dir().z() );
      shwvstrue->Fill( evt->Get_MC_energy(), evt->Get_shower_energy() );
      shwres->Fill( evt->Get_MC_energy(), ( evt->Get_MC_energy() - evt->Get_shower_energy() )/evt->Get_MC_energy() * 100 );
    }

  }

  // create projections
  vector<TH1D*> ptrk, pshw;
  
  for (Int_t xbin = 1; xbin <= trkres->GetXaxis()->GetNbins(); xbin++) {

    TString name = "E=" + (TString)to_string( trkres->GetXaxis()->GetBinCenter(xbin) ) + "GeV";
    
    ptrk.push_back( trkres->ProjectionY("trk_" + name, xbin, xbin) );
    pshw.push_back( shwres->ProjectionY("shw_" + name, xbin, xbin) );

    ptrk.back()->SetTitle(name);
    
  }

  // normalise the e_vs_true hists
  vector<TH2D*> hs = {trkvstrue, shwvstrue};
  for (auto h: hs) {

    h->GetXaxis()->SetTitle("True energy [GeV]");

    for (Int_t xbin = 1; xbin <= h->GetXaxis()->GetNbins(); xbin++) {

      TH1D* proj = h->ProjectionY("proj",xbin,xbin);
      Double_t I = proj->Integral();
      delete proj;

      for (Int_t ybin = 1; ybin <= h->GetYaxis()->GetNbins(); ybin++) {
	Double_t bc = h->GetBinContent(xbin,ybin);
	h->SetBinContent(xbin, ybin, bc/I);
      }

    }

  }
  
  
  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->DivideSquare( ptrk.size() );

  for (Int_t i = 0; i < ptrk.size(); i++) {
    c1->cd(i+1);
    ptrk[i]->Draw();
    pshw[i]->SetLineColor(kRed);
    pshw[i]->Draw("same");
  }

  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->Divide(2,2);
  c2->cd(1);
  trkvstrue->GetYaxis()->SetTitle("track reco energy");
  trkvstrue->Draw("colz");
  c2->cd(2);
  shwvstrue->GetYaxis()->SetTitle("shower reco energy");
  shwvstrue->Draw("colz");
  c2->cd(3);
  trkdir->Draw("colz");
  c2->cd(4);
  shwdir->Draw("colz");
  
}
