#include "EffMass.h"
#include "EMextr.h"
#include "DetResponse.h"
#include "NMHUtils.h"
#include <iostream>
#include "TRandom.h"
#include "TMath.h"
#include "TFile.h"
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

using namespace std;

/** 
    This applications tests the extrapolation in effective mass.

    \return 0 if test worked, 1 otherwise.
*/
int main(const int argc, const char** argv) {

  TString outfile;

  try {
    JParser<> zap("This application tests the effective mass extrapolator class `EXextr` against the purely MC statistics based effective mass class `EffMass`");
    zap['o'] = make_field(outfile, "Draw canvases for visual inspection of the effective mass curves into this file") = "";
    if ( zap.read(argc, argv) != 0 ) return 1;
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  Int_t ret = 0;

  //=============================================================================
  // fetch an effective mass file from the server
  //=============================================================================

  TString effmfile = "EffMass_ORCA115_20x9m_ECAP190222.root";
  TString syscmd   = "wget http://sftp.km3net.de/MONA/" + effmfile;
  Int_t sysret     = system( syscmd );

  if (sysret != 0) {
    cout << "ERROR! effmass: file retrieval returned " << sysret << ", exiting" << endl;
    system("rm " + effmfile);
    return 1;
  }
  else {
    cout << "NOTICE effmass: file retrieval done" << endl;
  }

  //=============================================================================
  // fetch an effective mass file from the server
  //=============================================================================

  // response where true and reco have same binning, should yield identical results
  DetResponse R1(DetResponse::track, "R1", 40, 1, 100, 40, -1, 1, 2, 0, 1);

  // init the standard effective mass calculator
  EffMass em(effmfile, 40, 40, 2);

  // init the interpolated eff mass calculator with same binning, results should match identically
  EMextr  ex1(R1.GetHist3DTrue(), effmfile);
  
  TRandom3 rand(0);

  //==========================================================================
  // test 1 - see that em and ex1 give identical results, as same binning is used
  //==========================================================================

  Int_t trials = 1000;

  for (Int_t N = 0; N < trials; N++) {

    Double_t f    = rand.Integer(3);
    Double_t iscc = rand.Integer(2);
    Double_t isnb = rand.Integer(2);
    Double_t E    = rand.Uniform(1.0, 100.0);
    Double_t ct   = rand.Uniform(-1.0, 1.0);
    Double_t by   = rand.Uniform(0.0, 1.0);

    Double_t effm1 = em.GetMeff(f, iscc, isnb, E, ct, by);
    Double_t effm2 = ex1.GetMeff(f, iscc, isnb, E, ct, by);

    Double_t diff = TMath::Abs(effm1 - effm2)/effm1 * 100;
    if ( diff > 1e-3 ) {
      cout << "ERROR! effmass: effective masses differ, values and diff: " << effm1 << "\t" << effm2 << "\t" << diff << "%" << endl;
      ret += 1;
    }

  }

  //==========================================================================
  // init a response with different binning and compare interpolated results
  // otherwise I would be comparing results in different bins, which are expected to differ
  //==========================================================================

  // desponse where true and reco have different binning and true energy extends MC range
  DetResponse R2(DetResponse::track, "R2", 70, 0.1, 150, 70, -1, 1, 2, 0, 1, 40, 1, 100, 40, -1, 1, 1, 0, 1);

  // init the interpolated eff mass calculator with different binning, results should match within tolerance
  EMextr  ex2(R2.GetHist3DTrue(), effmfile);

  Double_t ave_diff = 0.0;

  for (Int_t N = 0; N < trials; N++) {

    Double_t f    = rand.Integer(3);
    Double_t iscc = rand.Integer(2);
    Double_t isnb = rand.Integer(2);

    // only take data points in MC interpolation range
    Double_t E    = rand.Uniform(10.0, 50.0);
    Double_t ct   = rand.Uniform(-0.9, 0.9);
    Double_t by   = rand.Uniform(0.0, 1.0);

    Double_t effm1 = em.GetMeff(f, iscc, isnb, E, ct, by, kTRUE);

    // for EMextr I need to do the interpolation step here manually
    TH3D* h3d = (TH3D*)ex2.GetMeff3DH(f, iscc, isnb)->Clone("tmp3d");
    Int_t byb = h3d->GetZaxis()->FindBin(by);
    h3d->GetZaxis()->SetRange(byb, byb);
    TH2D* h2d = (TH2D*)h3d->Project3D("yx");
    
    Double_t effm2 = h2d->Interpolate(E, ct);

    delete h2d;
    delete h3d;

    Double_t diff = TMath::Abs(effm1 - effm2)/effm1 * 100;
    ave_diff += diff;

  }

  if ( ave_diff/trials > 3.0 ) {
    cout << "NOTICE effmass: average difference " << ave_diff/trials << " is larger than 3.0%" << endl;
    ret += 1;
  }
  
  //==========================================================================
  // draw plots cc-nu plots for each flavor, optinally
  //==========================================================================

  if ( outfile != "" ) {

    TFile fout(outfile, "RECREATE");

    vector<Int_t> flavs = {0, 1, 2};

    for (auto f: flavs) {

      vector<TH1D*> slices_em, slices_extr;

      for (Double_t ct = -0.95; ct <= 0.95; ct += 0.2) {
      
	TString nt_em   = "slice_em_flavor_" + (TString)to_string(f) + "_ct_" + (TString)to_string(ct);
	TString nt_extr = "slice_extr_flavor_" + (TString)to_string(f) + "_ct_" + (TString)to_string(ct);

	TH1D* slice_em = em.GetSlice(1, 1, 0, ct, 0.25);
	slice_em->SetNameTitle( nt_em, nt_em );

	TH1D* slice_extr = ex2.GetSlice(1, 1, 0, ct, 0.25);
	slice_extr->SetNameTitle( nt_extr, nt_extr );
      
	slice_em->SetLineColor(kRed);
	slice_em->SetLineWidth(2);

	slice_extr->SetLineColor(kBlue);
	slice_extr->SetLineWidth(2);

	slices_em.push_back( slice_em );
	slices_extr.push_back( slice_extr );
	
      }

      TString cname = "cflav_" + (TString)to_string( f );
      TCanvas *c = new TCanvas(cname, cname, 1);
      c->DivideSquare( slices_em.size() );
      for (Int_t i = 0; i < (Int_t)slices_em.size(); i++) {
	c->cd(i+1);
	slices_extr[i]->Draw("HIST");
	slices_em[i]->Draw("HISTsame");
      }
      c->Write();

    }

    fout.Close();

  }

  //==========================================================================
  // clean-up
  //==========================================================================
  
  sysret = system("rm EffMass_ORCA115_20x9m_ECAP190222.root*");
  if ( sysret != 0 ) {
    cout << "NOTICE effmass: file deletion returned " << sysret << endl;
  }

  cout << "NOTICE: effmass finished" << endl;
  return ret;

}
