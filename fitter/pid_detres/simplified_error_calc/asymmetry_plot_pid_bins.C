/* This plot is a graph with error bars of the asymmetry for track, shower and combined
 * as a function of PID number of bins, ranging from 2 to 40.
 *
 * The version for simplified errors is in ./simplified_error_calc/
 */

void asymmetry_plot_pid_bins() {
  TCanvas *c1 = new TCanvas("c1", "c1", 200, 10, 700, 500);

  c1->SetGrid();
  c1->GetFrame()->SetBorderSize(12);

  const Int_t n = 9;
  // combined
  Float_t x_c[n]  = {2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0};
  Float_t y_c[n]  = {6.09673, 6.79003, 7.4376, 7.83659, 8.34702, 8.36275, 8.96836, 9.39361, 9.80033};
  Float_t ex_c[n] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Float_t ey_c[n] = {0.0149421, 0.0181562, 0.244886, 0.233527, 0.306792, 0.221089, 0.207093, 0.274468, 0.323194};
  TGraphErrors *gc = new TGraphErrors(n, x_c, y_c, ex_c, ey_c);
  gc->Draw("AP");
  gc->SetTitle("Asymmetry as function of number of PID bins");
  gc->GetXaxis()->SetTitle("Number of bins used");
  gc->GetYaxis()->SetTitle("Asymmetry");
  gc->SetMarkerColor(kYellow+2);
  gc->SetMarkerStyle(4);
  gc->SetMarkerSize(1);
  gc->GetYaxis()->SetRangeUser(0,13);
  c1->Update();

  // track
  Float_t x_t[n]  = {2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0};
  Float_t y_t[n]  = {3.43853, 3.79568, 4.47484, 4.79428, 5.11288, 5.00107, 5.69825, 5.9665, 6.33294};
  Float_t ex_t[n] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Float_t ey_t[n] = {0.0232898, 0.0284127, 0.406697, 0.379761, 0.35851, 0.367816, 0.324082, 0.311018, 0.412574};
  TGraphErrors *gt = new TGraphErrors(n, x_t, y_t, ex_t, ey_t);
  gt->Draw("P");
  gt->SetMarkerColor(kBlue+2);
  gt->SetMarkerStyle(4);
  c1->Update();

  // shower
  Float_t x_s[n]  = {2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0};
  Float_t y_s[n]  = {5.03454, 5.63004, 5.94086, 6.19895, 6.59782, 6.7026, 6.92541, 7.2554, 7.47933};
  Float_t ex_s[n] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  Float_t ey_s[n] = {0.00862511, 0.0106089, 0.0122967, 0.0298512, 0.271033, 0.0278348, 0.0285863, 0.246701, 0.239386};
  TGraphErrors *gs = new TGraphErrors(n, x_s, y_s, ex_s, ey_s);
  gs->Draw("P");
  gs->SetMarkerColor(kRed+2);
  gs->SetMarkerStyle(4);
  c1->Update();

  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->SetHeader("Channel","C"); // option "C" allows to center the header
  legend->AddEntry(gt, "Track", "ep");
  legend->AddEntry(gs, "Shower", "ep");
  legend->AddEntry(gc, "Combined", "ep");
  legend->Draw();

}
