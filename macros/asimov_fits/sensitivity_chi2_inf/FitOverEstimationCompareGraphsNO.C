#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"

void FitOverEstimationCompareGraphsNO() {

  std::vector<Int_t> pid_cats = {2,3,4,5,10};

  TFile *f = TFile::Open("data_chi2_th23_extraplotation.root", "READ");
    
  TString s_tree; // tree name
  TString plot_selection; // branches to use in pot
  TString plot_title;
  TString legend_label;
  TString legend_plotstyle;
  TString out_name;

  Bool_t plotInf[2] = {kFALSE, kTRUE};

  for (auto p: plotInf) {
    Bool_t plot_inf = p;

    TCanvas* c1 = new TCanvas("c1", "c1");
    TLegend* leg = new TLegend(0.1, 0.7, 0.4, 0.9);

    if (plot_inf) {
      s_tree = "inf_";
      plot_selection = "th23:sqrt_inf_chi2:sqrt_inf_err";
      plot_title = "#Delta #chi^{2} extrapolated to infinite statistics";
      legend_label = "NO, %i PID categories";
      legend_plotstyle = "lpe";
      out_name = "inf";
    } else {
      s_tree = "chi2_";
      plot_selection = "th23:sqrt_fit_chi2";
      plot_title = "#Delta #chi^{2} at finite statistics";
      legend_label = "NO, %i PID categories";
      legend_plotstyle = "lp";
      out_name = "fin";
    }

    TTree* t_2  = (TTree*)f->Get(s_tree + "2bins_io");
    TTree* t_3  = (TTree*)f->Get(s_tree + "3bins_io");
    TTree* t_4  = (TTree*)f->Get(s_tree + "4bins_io");
    TTree* t_5  = (TTree*)f->Get(s_tree + "5bins_io");
    TTree* t_10 = (TTree*)f->Get(s_tree + "10bins_io");
  
    // Draw the fitted chi2_inf
  
    Int_t n2 = t_2->Draw(plot_selection, "", "goff");
    TGraphErrors *g2 = new TGraphErrors(n2, t_2->GetV1(), t_2->GetV2(), 0, t_2->GetV3()); 
    g2->Sort(&TGraph::CompareX);
    Int_t n3 = t_3->Draw(plot_selection, "", "goff");
    TGraphErrors *g3 = new TGraphErrors(n3, t_3->GetV1(), t_3->GetV2(), 0, t_2->GetV3()); 
    g3->Sort(&TGraph::CompareX);
    Int_t n4 = t_4->Draw(plot_selection, "", "goff");
    TGraphErrors *g4 = new TGraphErrors(n4, t_4->GetV1(), t_4->GetV2(), 0, t_2->GetV3()); 
    g4->Sort(&TGraph::CompareX);
    Int_t n5 = t_5->Draw(plot_selection, "", "goff");
    TGraphErrors *g5 = new TGraphErrors(n5, t_5->GetV1(), t_5->GetV2(), 0, t_2->GetV3()); 
    g5->Sort(&TGraph::CompareX);
    //Int_t n10 = t_10->Draw(plot_selection, "", "goff");
    //TGraphErrors *g10 = new TGraphErrors(n10, t_10->GetV1(), t_10->GetV2(), 0, t_2->GetV3()); 
    //g10->Sort(&TGraph::CompareX);
  
    g2->SetLineColor(kBlue+1);
    g2->SetMarkerColor(kBlue+1);
  //  g2->SetMarkerStyle(21);
  
    g3->SetLineColor(kGreen+3);
    g3->SetMarkerColor(kGreen+3);
  //  g3->SetMarkerStyle(21);
  
    g4->SetLineColor(kRed+2);
    g4->SetMarkerColor(kRed+2);
  //  g4->SetMarkerStyle(21);
  
    g5->SetLineColor(kAzure+1);
    g5->SetMarkerColor(kAzure+1);
  //  g5->SetMarkerStyle(21);
  
    //g10->SetLineColor(kBlack);
    //g10->SetMarkerColor(kBlack);
  //  g10->SetMarkerStyle(21);
  
    g2->SetMinimum(0);
    g2->SetMaximum(7);
    g2->Draw("apl"); 
    g3->Draw("pl"); 
    g4->Draw("pl"); 
    g5->Draw("pl"); 
    //g10->Draw("pl"); 
  
    g2->SetTitle(plot_title);
    g2->GetXaxis()->SetTitle("#theta_{23}");
    g2->GetYaxis()->SetTitle("#sqrt{ #Delta #chi^{2} }");
    g2->GetXaxis()->SetTitleSize(0.04);
    g2->GetYaxis()->SetTitleSize(0.04);
  
    leg->AddEntry(g2, Form(legend_label, 2), legend_plotstyle);
    leg->AddEntry(g3, Form(legend_label, 3), legend_plotstyle);
    leg->AddEntry(g4, Form(legend_label, 4), legend_plotstyle);
    leg->AddEntry(g5, Form(legend_label, 5), legend_plotstyle);
    //leg->AddEntry(g10, Form(legend_label, 10), legend_plotstyle);
    leg->Draw();
  
    TFile* fout = new TFile("FitOverEstimationCompareGraphs.root", "UPDATE");
    fout->cd();
    g2->Write("comaprison_2bins_no_" + out_name);
    g3->Write("comaprison_3bins_no_" + out_name);
    g4->Write("comaprison_4bins_no_" + out_name);
    g5->Write("comaprison_5bins_no_" + out_name);
    fout->Close();

    if (plot_inf) c1->SaveAs("CompareChi2AtInfStatisticsNO.pdf");

    delete c1;
    delete leg;
    delete fout; // Clear memory for second iteration.
  }
}
