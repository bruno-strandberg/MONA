void summaryparser() {

  // first create some pseudodata that can be read in later

  TString fname = "dataparsing.root";

  SummaryParser wr(fname, kFALSE);
  for (Int_t i = 0; i < 100; i++) {
    wr.GetEvt()->FillPseudoData();
    wr.GetTree()->Fill();
  }
  wr.WriteAndClose();

  // now read in the data using summaryparser and cout some variables

  SummaryParser sp(fname);

  for (Int_t i = 0; i < sp.GetTree()->GetEntries(); i++) {
    cout << "MC Energy, dir_z in summary: " 
	 << sp.GetEvt(i)->Get_MC_energy() << "\t" 
	 << sp.GetEvt(i)->Get_MC_dir().z() << endl;
  }
  
  // finally remove the pseudodata file

  system("rm " + fname);
}
