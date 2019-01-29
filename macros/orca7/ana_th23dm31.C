using namespace RooFit;

/* 
   This macro analysis the output of `th23dm31`. It plots the true th23, dm31 point and the values extracted from fits to pseudo-experiments on a TGraph.

   infile  - output of th23dm31.C macro
   outfile - where output graphs should be written. By default empty string and outfile not created.
*/
void ana_th23dm31(TString infile, TString outfile="") {

  //----------------------------------------------------------------------------
  // read in the true value
  //----------------------------------------------------------------------------

  FileHeader head("ana_th23dm31");
  head.ReadHeader(infile);

  Double_t th23 = stod( (string)head.GetParameter("SinsqTh23") );
  Double_t dm31 = stod( (string)head.GetParameter("Dm31") );

  TGraph *trueplot = new TGraph();
  trueplot->SetPoint(0, th23, dm31);

  //----------------------------------------------------------------------------
  // loop over experiment results and add points to results plot
  //----------------------------------------------------------------------------

  TGraph *resplot = new TGraph();
  TH1D* hth23 = new TH1D("hth23","hth23", 50, 0.25, 0.75);
  TH1D* hdm31 = new TH1D("hdm31","hdm31", 30, 0.0015, 0.0045);
  
  TFile fin(infile,"READ");
  
  TObject *obj;
  TKey    *key;
  TIter    next( fin.GetListOfKeys() );

  while ( (key = (TKey *) next() ) ) {

    obj = fin.Get( key->GetName() );
    TString objname = obj->GetName();

    if ( objname.Contains("result_") ) {
      
      RooFitResult *res = (RooFitResult*)obj;

      Double_t th23fit = ( (RooRealVar*)res->floatParsFinal().find("SinsqTh23") )->getVal();
      Double_t dm31fit = ( (RooRealVar*)res->floatParsFinal().find("Dm31") )->getVal();

      resplot->SetPoint( resplot->GetN(), th23fit, dm31fit );
      hth23->Fill(th23fit);
      hdm31->Fill(dm31fit);

    }
  }

  TLine *trueth23 = new TLine( th23, 0, th23, hth23->GetMaximum() );
  trueth23->SetLineColor(kRed);
  trueth23->SetLineWidth(2);
  TLine *truedm31 = new TLine( dm31, 0, dm31, hdm31->GetMaximum() );
  truedm31->SetLineColor(kRed);
  truedm31->SetLineWidth(2);

  //----------------------------------------------------------------------------
  // plot stuff
  //----------------------------------------------------------------------------
  trueplot->SetMarkerStyle(20);
  trueplot->SetMarkerColor(kRed);
  trueplot->SetMarkerSize(2);
  
  resplot->SetMarkerStyle(5);
  resplot->SetMarkerColor(kBlue);
  
  TCanvas *c1 = new TCanvas("c1","c1",1);
  resplot->Draw("AP");
  trueplot->Draw("sameP");
  
  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->Divide(2,1);
  c2->cd(1);
  hth23->Draw();
  trueth23->Draw();
  c2->cd(2);
  hdm31->Draw();
  truedm31->Draw();  

  //----------------------------------------------------------------------------
  // write plots to output
  //----------------------------------------------------------------------------

  if (outfile != "") {
    TFile fout(outfile, "RECREATE");
    trueplot->Write("truepoint");
    resplot->Write("pseudoexps");
    hdm31->Write();
    hth23->Write();
    fout.Close();
  }
 
}
