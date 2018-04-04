#include "GSGParser.h"
#include<iostream>
using namespace std;

void gsgparser() {

  gSystem->Load("$NMHDIR/common_software/libnmhsoft.so");

  GSGParser gp1("../data/mc_start/data_atmnu/gSeaGen_muon-CC_3-100GeV_1.root");
  GSGParser gp2("../data/testing/gSeaGen_muon-CC_3-100GeV_1.evt");

  cout << "********************HEADER DATA****************************" << endl;
  cout << "root file\tevt file" << endl;
  cout << gp1.fNgen << "\t" << gp2.fNgen << endl;
  cout << gp1.fVint << "\t" << gp2.fVint << endl;
  cout << gp1.fRho_seawater << "\t" << gp2.fRho_seawater << endl;
  cout << gp1.fVcan << "\t" << gp2.fVcan << endl;
  cout << gp1.fRcan << "\t" << gp2.fRcan << endl;
  cout << gp1.fZcan_min << "\t" << gp2.fZcan_min << endl;
  cout << gp1.fZcan_max << "\t" << gp2.fZcan_max << endl;
  cout << gp1.fAgen << "\t" << gp2.fAgen << endl;
  cout << gp1.fE_min << "\t" << gp2.fE_min << endl;
  cout << gp1.fE_max << "\t" << gp2.fE_max << endl;
  cout << gp1.fCt_min << "\t" << gp2.fCt_min << endl;
  cout << gp1.fCt_max << "\t" << gp2.fCt_max << endl;
  cout << gp1.fE_power << "\t" << gp2.fE_power << endl;


  cout << "********************TEST DIFFERENT LOOP METHODS************" << endl;
  
  TH1D *h1 = new TH1D("h1","h1",100,0,100);
  TH1D *h2 = new TH1D("h2","h2",100,0,100);

  while ( gp1.NextEvent() ) {
    h1->Fill( gp1.Neutrino_E );
  }

  for (Int_t i = 0; i < gp1.fChain->GetEntries(); i++) {
    gp1.fChain->GetEntry(i);
    h2->Fill( gp1.Neutrino_E );
  }
  
  cout << "********************TEST DIFFERENT FILE EXTENSIONS*********" << endl;

  TH1D *h3 = new TH1D("h3","h3",100,0,100);
  TH1D *h4 = new TH1D("h4","h4",100,0,100);

  Int_t i = 0;
  while (i < 3) {
    gp1.fChain->GetEntry(i);
    gp2.NextEvent();
    i++;
    cout << gp1.Neutrino_V1 << "\t" << gp2.Neutrino_V1 << endl;
    cout << gp1.Neutrino_V2 << "\t" << gp2.Neutrino_V2 << endl;
    cout << gp1.Neutrino_V3 << "\t" << gp2.Neutrino_V3 << endl;
    cout << gp1.Neutrino_D1 << "\t" << gp2.Neutrino_D1 << endl;
    cout << gp1.Neutrino_D2 << "\t" << gp2.Neutrino_D2 << endl;
    cout << gp1.Neutrino_D3 << "\t" << gp2.Neutrino_D3 << endl;
    cout << gp1.Neutrino_E << "\t" << gp2.Neutrino_E << endl;
    cout << gp1.Neutrino_PdgCode << "\t" << gp2.Neutrino_PdgCode << endl;
    cout << "=========================================================" << endl;
  }

  for (Int_t j = i; j < gp1.fChain->GetEntries(); j++) {
    gp1.GetEntry(j);
    h3->Fill(gp1.Neutrino_E);
  }

  while ( gp2.NextEvent() ) {
    h4->Fill( gp2.Neutrino_E );
  }

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->Divide(2,2);
  c1->cd(1);
  h1->Draw();
  c1->cd(2);
  h2->Draw();
  c1->cd(3);
  h3->Draw();
  c1->cd(4);
  h4->Draw();

}
