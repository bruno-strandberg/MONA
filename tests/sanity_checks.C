void sanity_checks(TString fname) {

  TFile *f = new TFile(fname, "READ");
  TTree *t = (TTree*)f->Get("PID");

  //this should have exactly 0 entries
  TH1D *h1 = new TH1D("h1","h1",200, -1, 1);
  t->Draw("gandalf_dir_z>>h1","gandalf_dir_z<0&&dusj_dir_z<0&&recolns_dir_z<0");

  //these histograms should have at least one quantile empty
  TH2D *h2 = new TH2D("h2","h2",200,-1,1,200,-1,1);
  TH2D *h3 = new TH2D("h3","h3",200,-1,1,200,-1,1);
  TH2D *h4 = new TH2D("h4","h4",200,-1,1,200,-1,1);
  t->Draw("gandalf_dir_z:dusj_dir_z>>h2","recolns_dir_z<0");
  t->Draw("gandalf_dir_z:recolns_dir_z>>h3","dusj_dir_z<0");
  t->Draw("dusj_dir_z:recolns_dir_z>>h4","gandalf_dir_z<0");

  //these histograms I thought should have only dir_z > 0, but that is not the case
  TH1D *h5 = new TH1D("h5","h5",200,-1,1);
  TH1D *h6 = new TH1D("h6","h6",200,-1,1);
  TH1D *h7 = new TH1D("h7","h7",200,-1,1);
  t->Draw("gandalf_dir_z>>h5","gandalf_is_good>0");
  t->Draw("dusj_dir_z>>h6","dusj_is_good>0");
  t->Draw("recolns_dir_z>>h7","recolns_is_good>0");

  //these histograms I should have only dir_z > 0, as of 14.03.2018
  TH1D *h8  = new TH1D("h8","h8",200,-1,1);
  TH1D *h9  = new TH1D("h9","h9",200,-1,1);
  TH1D *h10 = new TH1D("h10","h10",200,-1,1);
  t->Draw("gandalf_dir_z>>h8","gandalf_is_good>0&&gandalf_is_selected>0");
  t->Draw("dusj_dir_z>>h9","dusj_is_good>0&&dusj_is_selected>0");
  t->Draw("recolns_dir_z>>h10","recolns_is_good>0&&recolns_is_selected>0");
  
  TFile *fout = new TFile("~/Dropbox/Working/sanity_out.root","RECREATE");
  h1->Write();
  h2->Write();
  h3->Write();
  h4->Write();
  h5->Write();
  h6->Write();
  h7->Write();
  h8->Write();
  h9->Write();
  h10->Write();
  fout->Close();
  
}
