#include "FileHeader.h"
#include "TH2.h"
#include "TFile.h"
#include "TGraph.h"

/**
 *  This macro uses the output of EffMhists.C to produce effective mass plots.
 *
 * \param  effmhists_file  Output of EffMhists.C macro.
 * \param  outname         Name of the file where outputs are written
 * \param  rebinX          Rebin X (energy) axis.
 * \param  rebinY          Rebin Y (costheta) axis.
 * \param  upgoing         For angle-averaged Meff curves only consider upgoing
 *
 */

void EffMass(TString effmhists_file, TString outname, 
	     Int_t rebinX = 3, Int_t rebinY = 5, Bool_t upgoing = kTRUE) {

  //-----------------------------------------
  //Get the histograms from input, rebin, divide, rename
  //-----------------------------------------
  FileHeader emh("emhheader");                            //read the header from effmhists output
  emh.ReadHeader(effmhists_file);
  
  FileHeader h("EffMass");                                //create header for this application
  h.AddParameter( "effmhists_file", effmhists_file );     //add the input and output names
  h.AddParameter( "outname", outname );

  // manually add some header fields from EffMhists header 
  // if effmhists_file is hadd, there are multiple entries for these - adding them manually
  // makes the header of this application shorter
  h.AddParameter( emh, "Rvol" );
  h.AddParameter( emh, "Zmin" );
  h.AddParameter( emh, "Zmax" );
  h.AddParameter( emh, "atmmu_cut" );
  h.AddParameter( emh, "noise_cut" );

  TFile *fin  = new TFile(effmhists_file, "READ");

  TH2D *h_gen_nu          = (TH2D*)fin->Get("Generated_scaled_nu");
  TH2D *h_gen_nub         = (TH2D*)fin->Get("Generated_scaled_nub");
  TH2D *h_det_nu          = (TH2D*)fin->Get("Detected_nu");
  TH2D *h_det_nub         = (TH2D*)fin->Get("Detected_nub");

  h_gen_nu         ->Rebin2D(rebinX, rebinY);
  h_gen_nub        ->Rebin2D(rebinX, rebinY);
  h_det_nu         ->Rebin2D(rebinX, rebinY);          
  h_det_nub        ->Rebin2D(rebinX, rebinY);         

  h_det_nu         ->Divide(h_gen_nu);
  h_det_nub        ->Divide(h_gen_nub);
  
  h_det_nu         ->SetNameTitle("Meff_nu","Meff_nu");
  h_det_nub        ->SetNameTitle("Meff_nub","Meff_nub");

  //-----------------------------------------
  //create Meff curves
  //-----------------------------------------
  TGraph *g_nu  = new TGraph();
  TGraph *g_nub = new TGraph();

  Int_t ct_max_bin = h_det_nu->GetYaxis()->GetNbins();
  if (upgoing) ct_max_bin = h_det_nu->GetYaxis()->FindBin(-0.0001);
  
  for (Int_t ebin = 1; ebin <= h_det_nu->GetXaxis()->GetNbins(); ebin++) {

    Double_t int_nu          = 0;
    Double_t int_nub         = 0;
    Double_t nbins           = 0;

    for (Int_t ctbin = 1; ctbin <= ct_max_bin; ctbin++) {
      int_nu          += h_det_nu ->GetBinContent(ebin, ctbin);
      int_nub         += h_det_nub->GetBinContent(ebin, ctbin);
      nbins++;
    }

    Int_t point     = g_nu->GetN();
    Double_t energy = h_det_nu->GetXaxis()->GetBinCenter(ebin);
    g_nu ->SetPoint(point, energy, int_nu/nbins  );
    g_nub->SetPoint(point, energy, int_nub/nbins );

  }

  TString title = "cos#theta average Meff [Tons]";
  if (upgoing) title += ", upgoing";
  g_nu ->SetNameTitle("E_vs_Meff_nu" ,title);
  g_nub->SetNameTitle("E_vs_Meff_nub",title);
  g_nu->SetLineColor(kRed);
  g_nub->SetLineColor(kBlue);
  
  //-----------------------------------------
  //write out
  //-----------------------------------------
  TFile *fout = new TFile(outname, "RECREATE");
  h_det_nu->Write();
  h_det_nub->Write();  
  g_nu->Write();
  g_nub->Write();
  h.WriteHeader(fout);
  fout->Close();

  //-----------------------------------------
  //cleanup
  //-----------------------------------------
  fin->Close();

  if (fin)   delete fin;
  if (fout)  delete fout;
  if (g_nu)  delete g_nu;
  if (g_nub) delete g_nub;

}
