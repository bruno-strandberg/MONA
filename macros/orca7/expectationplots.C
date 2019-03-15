#include "NMHUtils.h"


/* 
   In this script ORCA 7-line events are distributed to tracks and showers, using track score value 0.6, as in ORCA 115 production from 2015. This script creates:
   
   1. plots that depict the number of detected neutrinos, muons and noise in 1 year with ORCA 7-line detector in energy dimension.
   2. plots that depict the number of detected neutrinos in E-vs-costh plot, for tracks and showers.
   3. plots that depict the "asymmetry" in the track channel and in the shower channel for two different value pairs of (th23, dm31).

*/

void expectationplots() {

  //------------------------------------------------------------
  // input parameters
  //------------------------------------------------------------

  Int_t    rt    =   1;    //1 year
  Double_t emin  =   1;    // GeV min
  Double_t emax  = 100;    // GeV max
  Double_t ctmin =  -1;
  Double_t ctmax =   1;
  Double_t bymin =   0;
  Double_t bymax =   1;
  TString dataf  = (TString)getenv("MONADIR") + "/data/ORCA_MC_summary_ORCA7_23x9m_ECAP1018.root";
  TString effmf  = (TString)getenv("MONADIR") + "/data/eff_mass/EffMass_ORCA7_23x9m_ECAP1018.root";

  Double_t muoncut  = 0.05;
  Double_t noisecut = 0.01;
  Double_t trackcut = 0.6;

  //------------------------------------------------------------
  // initialise and fill track and shower responses
  //------------------------------------------------------------

  DetResponse track(DetResponse::track, "track");
  track.AddCut( &SummaryEvent::Get_track_ql0       , std::greater<double>()   ,      0.5, true);
  track.AddCut( &SummaryEvent::Get_RDF_noise_score , std::less_equal<double>(), noisecut, true);
  track.AddCut( &SummaryEvent::Get_RDF_muon_score  , std::less_equal<double>(),  muoncut, true);
  track.AddCut( &SummaryEvent::Get_RDF_track_score , std::greater<double>()   , trackcut, true);

  DetResponse shower(DetResponse::shower, "shower");
  shower.AddCut( &SummaryEvent::Get_shower_ql0     , std::greater<double>()   ,      0.5, true);
  shower.AddCut( &SummaryEvent::Get_RDF_noise_score, std::less_equal<double>(), noisecut, true);
  shower.AddCut( &SummaryEvent::Get_RDF_muon_score , std::less_equal<double>(),  muoncut, true);
  shower.AddCut( &SummaryEvent::Get_RDF_track_score, std::less_equal<double>(), trackcut, true);

  SummaryParser sp(dataf);
  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    if (i % 100000 == 0) cout << "Event: " << i << " out of: " << sp.GetTree()->GetEntries() << endl;
    track.Fill ( sp.GetEvt(i) );
    shower.Fill( sp.GetEvt(i) );
  }

  //------------------------------------------------------------
  // initialise the pdf's
  //------------------------------------------------------------

  FitUtil futil(rt, track.GetHist3D(), emin, emax, ctmin, ctmax, bymin, bymax, effmf);
  FitPDF trackpdf ("trackpdf" ,"trackpdf" , &futil, &track);
  FitPDF showerpdf("showerpdf","showerpdf", &futil, &shower);

  //------------------------------------------------------------
  // plot the expected number of neutrinos, muons and noise 
  //------------------------------------------------------------
  gStyle->SetOptStat(0);
  vector<FitPDF*> pdfs = {&trackpdf, &showerpdf};

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->Divide(2,1);
  Int_t pad = 1;

  for (auto pdf: pdfs) {

    TH1D* nus   = (TH1D*)pdf->GetExpValHist()->Project3D("x");
    TH1D* muons = (TH1D*)pdf->GetResponse()->GetHistAtmMu1y()->Project3D("x");
    TH1D* noise = (TH1D*)pdf->GetResponse()->GetHistNoise1y()->Project3D("x");

    nus->SetLineColor(kRed);
    muons->SetLineColor(kBlack);
    noise->SetLineColor(kGreen);

    nus->GetXaxis()->SetTitle("Reco Energy [GeV]");
    nus->GetYaxis()->SetTitle("Events [1 year]");

    TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg->AddEntry(nus, "neutrinos", "l");
    leg->AddEntry(muons, "muons", "l");
    leg->AddEntry(noise, "noise", "l");

    c1->cd(pad);
    nus->Draw("HIST");
    noise->Draw("HISTsame");
    muons->Draw("HISTsame");
    leg->Draw();
    pad++;

    cout << "Neutrinos, muons and noise in 1 y for " << pdf->GetResponse()->Get_RespName() << ": " 
	 << nus->Integral() << "\t" << muons->Integral() << "\t" << noise->Integral() << endl;

  }
 
  //------------------------------------------------------------
  // create an asymmetry plot
  //------------------------------------------------------------
  RooRandom::randomGenerator()->SetSeed(416); // this seed controls the randomisation of osc parameters

  futil.GetVar("SinsqTh23")->randomize();
  futil.GetVar("Dm31")->randomize();

  cout << "*********************************************************************" << endl;
  cout << "1st pair: " << futil.GetVar("SinsqTh23")->getVal() << "\t" << futil.GetVar("Dm31")->getVal() << endl;
  cout << "*********************************************************************" << endl;

  TH2D* trks1 = (TH2D*)trackpdf.GetExpValHist()->Project3D("yx")->Clone();
  TH2D* shws1 = (TH2D*)showerpdf.GetExpValHist()->Project3D("yx")->Clone();

  futil.GetVar("SinsqTh23")->randomize();
  futil.GetVar("Dm31")->randomize();

  cout << "*********************************************************************" << endl;
  cout << "2nd pair: " << futil.GetVar("SinsqTh23")->getVal() << "\t" << futil.GetVar("Dm31")->getVal() << endl;
  cout << "*********************************************************************" << endl;

  TH2D* trks2 = (TH2D*)trackpdf.GetExpValHist()->Project3D("yx")->Clone();
  TH2D* shws2 = (TH2D*)showerpdf.GetExpValHist()->Project3D("yx")->Clone();
  
  // this is necessary because otherwise ROOT get's confused and the asymmetry calculator somehow
  // assumes the histograms to be the same...
  vector<TH2D*> hists = { trks1, shws1, trks2, shws2 }; 
  for (auto h: hists) h->SetDirectory(0);

  auto trk_asym = NMHUtils::Asymmetry(trks1, trks2, "trkasym");
  auto shw_asym = NMHUtils::Asymmetry(shws1, shws2, "shwasym");

  hists = { trks1, shws1, std::get<0>(trk_asym), std::get<0>(shw_asym) };
  
  for (auto h: hists) {
    h->GetXaxis()->SetTitle("Reco Energy [GeV]");
    h->GetYaxis()->SetTitle("Reco cos#theta");
    h->GetYaxis()->SetRangeUser(-1,0);
  }
  
  trks1->SetTitle("tracks, 1year");
  shws1->SetTitle("showers, 1year");
  std::get<0>(trk_asym)->SetTitle("tracks sensitivity, 1year");
  std::get<0>(shw_asym)->SetTitle("showers sensitivity, 1year");

  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->DivideSquare(4);
  c2->cd(1);
  trks1->Draw("colz");
  c2->cd(2);
  shws1->Draw("colz");
  c2->cd(3);
  std::get<0>(trk_asym)->Draw("colz");
  c2->cd(4);
  std::get<0>(shw_asym)->Draw("colz");

}
