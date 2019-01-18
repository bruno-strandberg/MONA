#include "NMHUtils.h"


/* Create 3D histograms with expectation values for 1 year of running for tracks and showers. I am performing some optimization of the muon and noise suppression cuts*/

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
  TString dataf  = "../../data/ORCA_MC_summary_ORCA7_23x9m_ECAP1018.root";
  TString effmf  = "../../data/eff_mass/EffMass_ORCA7_23x9m_ECAP1018.root";

  Double_t muoncut  = 0.05;
  Double_t noisecut = 0.01;
  Double_t trackcut = 0.6;

  //------------------------------------------------------------
  // initialise and fill track and shower responses; note all MC events are filled to both responses
  // at this time, i.e. I am not using track<->shower separation
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

    c1->cd(pad);
    nus->Draw("HIST");
    noise->Draw("HISTsame");
    muons->Draw("HISTsame");
    pad++;

    cout << "Neutrinos, muons and noise in 1 y for " << pdf->GetResponse()->Get_RespName() << ": " 
	 << nus->Integral() << "\t" << muons->Integral() << "\t" << noise->Integral() << endl;

  }
 
  //------------------------------------------------------------
  // create an asymmetry plot
  //------------------------------------------------------------
  TH2D* NOtrks = (TH2D*)trackpdf.GetExpValHist()->Project3D("yx")->Clone();
  TH2D* NOshws = (TH2D*)showerpdf.GetExpValHist()->Project3D("yx")->Clone();
  NOtrks->SetDirectory(0);
  NOshws->SetDirectory(0);

  futil.SetIOlims();
  
  TH2D* IOtrks = (TH2D*)trackpdf.GetExpValHist()->Project3D("yx")->Clone();
  TH2D* IOshws = (TH2D*)showerpdf.GetExpValHist()->Project3D("yx")->Clone();
  IOtrks->SetDirectory(0);
  IOshws->SetDirectory(0);

  auto trk_asym = NMHUtils::Asymmetry(NOtrks, IOtrks, "trkasym");
  auto shw_asym = NMHUtils::Asymmetry(NOshws, IOshws, "shwasym");

  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->Divide(2,1);
  c2->cd(1);
  std::get<0>(trk_asym)->Draw("colz");
  c2->cd(2);
  std::get<0>(shw_asym)->Draw("colz");

  TCanvas *c3 = new TCanvas("c3","c3",1);
  c3->Divide(2,2);
  c3->cd(1);
  NOtrks->Draw("colz");
  c3->cd(2);
  IOtrks->Draw("colz");
  c3->cd(3);
  NOshws->Draw("colz");
  c3->cd(4);
  IOshws->Draw("colz");

}
