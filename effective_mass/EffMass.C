
/**
 *  This macro uses the output of EffMhists.C to produce effective mass plots.
 *
 * \param  effmhists_file  Output of EffMhists.C macro.
 * \param  rebinX          Rebin X (energy) axis.
 * \param  rebinY          Rebin Y (costheta) axis.
 * \param  upgoing         For angle-averaged Meff curves only consider upgoing
 *
 */

void EffMass(TString effmhists_file, TString outname, 
	     Int_t rebinX = 2, Int_t rebinY = 5, Bool_t upgoing = kTRUE) {

  //-----------------------------------------
  //Get the histograms from input, rebin, divide, rename
  //-----------------------------------------
  TFile *fin  = new TFile(effmhists_file, "READ");

  TH2D *h_gen_nu          = (TH2D*)fin->Get("Generated_scaled_nu");
  TH2D *h_gen_nub         = (TH2D*)fin->Get("Generated_scaled_nub");
  TH2D *h_det_nu          = (TH2D*)fin->Get("Detected_nu");
  TH2D *h_det_nub         = (TH2D*)fin->Get("Detected_nub");
  TH2D *h_det_gandalf_nu  = (TH2D*)fin->Get("Detected_gandalf_nu");
  TH2D *h_det_gandalf_nub = (TH2D*)fin->Get("Detected_gandalf_nub");
  TH2D *h_det_shower_nu   = (TH2D*)fin->Get("Detected_shower_nu");
  TH2D *h_det_shower_nub  = (TH2D*)fin->Get("Detected_shower_nub");

  h_gen_nu         ->Rebin2D(rebinX, rebinY);
  h_gen_nub        ->Rebin2D(rebinX, rebinY);
  h_det_nu         ->Rebin2D(rebinX, rebinY);          
  h_det_nub        ->Rebin2D(rebinX, rebinY);         
  h_det_gandalf_nu ->Rebin2D(rebinX, rebinY);  
  h_det_gandalf_nub->Rebin2D(rebinX, rebinY); 
  h_det_shower_nu  ->Rebin2D(rebinX, rebinY);     
  h_det_shower_nub ->Rebin2D(rebinX, rebinY);  

  h_det_nu         ->Divide(h_gen_nu);
  h_det_nub        ->Divide(h_gen_nub);
  h_det_gandalf_nu ->Divide(h_gen_nu);
  h_det_gandalf_nub->Divide(h_gen_nub);
  h_det_shower_nu  ->Divide(h_gen_nu);
  h_det_shower_nub ->Divide(h_gen_nub);
  
  h_det_nu         ->SetNameTitle("Meff_nu","Meff_nu");
  h_det_nub        ->SetNameTitle("Meff_nub","Meff_nub");
  h_det_gandalf_nu ->SetNameTitle("Meff_gandalf_nu","Meff_gandalf_nu");
  h_det_gandalf_nub->SetNameTitle("Meff_gandalf_nub","Meff_gandalf_nub");
  h_det_shower_nu  ->SetNameTitle("Meff_shower_nu","Meff_shower_nu");
  h_det_shower_nub ->SetNameTitle("Meff_shower_nub","Meff_shower_nub");

  //-----------------------------------------
  //create Meff curves
  //-----------------------------------------
  TGraph *g_nu  = new TGraph();
  TGraph *g_nub = new TGraph();
  TGraph *g_gandalf_nu  = new TGraph();
  TGraph *g_gandalf_nub = new TGraph();
  TGraph *g_shower_nu   = new TGraph();
  TGraph *g_shower_nub  = new TGraph();

  Int_t ct_max_bin = h_det_nu->GetYaxis()->GetNbins();
  if (upgoing) ct_max_bin = h_det_nu->GetYaxis()->FindBin(-0.0001);
  
  for (Int_t ebin = 1; ebin <= h_det_nu->GetXaxis()->GetNbins(); ebin++) {

    Double_t int_nu          = 0;
    Double_t int_nub         = 0;
    Double_t int_gandalf_nu  = 0;
    Double_t int_gandalf_nub = 0;
    Double_t int_shower_nu   = 0;
    Double_t int_shower_nub  = 0;
    Double_t nbins           = 0;

    for (Int_t ctbin = 1; ctbin <= ct_max_bin; ctbin++) {
      int_nu          += h_det_nu ->GetBinContent(ebin, ctbin);
      int_nub         += h_det_nub->GetBinContent(ebin, ctbin);
      int_gandalf_nu  += h_det_gandalf_nu ->GetBinContent(ebin, ctbin);
      int_gandalf_nub += h_det_gandalf_nub->GetBinContent(ebin, ctbin);
      int_shower_nu   += h_det_shower_nu ->GetBinContent(ebin, ctbin);
      int_shower_nub  += h_det_shower_nub->GetBinContent(ebin, ctbin);
      nbins++;
    }

    Int_t point     = g_nu->GetN();
    Double_t energy = h_det_nu->GetXaxis()->GetBinCenter(ebin);
    g_nu ->SetPoint(point, energy, int_nu/nbins  );
    g_nub->SetPoint(point, energy, int_nub/nbins );
    g_gandalf_nu ->SetPoint(point, energy, int_gandalf_nu/nbins  );
    g_gandalf_nub->SetPoint(point, energy, int_gandalf_nub/nbins );
    g_shower_nu ->SetPoint(point, energy, int_shower_nu/nbins  );
    g_shower_nub->SetPoint(point, energy, int_shower_nub/nbins );

  }

  TString title = "cos#theta average Meff [Tons]";
  if (upgoing) title += ", upgoing";
  g_nu ->SetNameTitle("E_vs_Meff_nu" ,title);
  g_nub->SetNameTitle("E_vs_Meff_nub",title);
  g_gandalf_nu ->SetNameTitle("E_vs_Meff_gandalf_nu" ,title);
  g_gandalf_nub->SetNameTitle("E_vs_Meff_gandalf_nub",title);
  g_shower_nu ->SetNameTitle("E_vs_Meff_shower_nu" ,title);
  g_shower_nub->SetNameTitle("E_vs_Meff_shower_nub",title);
  g_nu->SetLineColor(kRed);
  g_nub->SetLineColor(kBlue);
  g_gandalf_nu->SetLineColor(kRed);
  g_gandalf_nub->SetLineColor(kBlue);
  g_shower_nu->SetLineColor(kRed);
  g_shower_nub->SetLineColor(kBlue);
  
  //-----------------------------------------
  //write out
  //-----------------------------------------
  TFile *fout = new TFile(outname, "RECREATE");
  h_det_nu         ->Write();
  h_det_nub        ->Write();  
  h_det_gandalf_nu ->Write();  
  h_det_gandalf_nub->Write();  
  h_det_shower_nu  ->Write();
  h_det_shower_nub ->Write();
  g_nu->Write();
  g_nub->Write();
  g_gandalf_nu->Write();
  g_gandalf_nub->Write();
  g_shower_nu->Write();
  g_shower_nub->Write();
  fout->Close();

  //-----------------------------------------
  //cleanup
  //-----------------------------------------
  fin->Close();

  if (fin)   delete fin;
  if (fout)  delete fout;
  if (g_nu)  delete g_nu;
  if (g_nub) delete g_nub;
  if (g_gandalf_nu)  delete g_gandalf_nu;
  if (g_gandalf_nub) delete g_gandalf_nub;
  if (g_shower_nu)   delete g_shower_nu;
  if (g_shower_nub)  delete g_shower_nub;

}
