#include "../common_software/SummaryParser.h"
#include "../common_software/GSGParser.h"
#include "../common_software/CScalculator.h"
#include "TSystem.h"
#include "TH2.h"
#include "TVector3.h"
#include "TMath.h"

//*********************************************************************
//functions

TString  ParseInputs(Int_t evt_type, Int_t int_type, Int_t en_low, Int_t run_nr);
void     InitHists();
Bool_t   VertexInCan(Double_t vx, Double_t vy, Double_t vz, Double_t Rcan, Double_t zcan_min, Double_t zcan_max);
void     FillDetected(Int_t veff_option);
void     FillGenerated_MCevts(Int_t veff_option);
void     FillGenerated_Agen(Int_t evt_type, Int_t int_type);

//*********************************************************************
//global histograms and variables

enum effvol_options { interaction_vol = 0, can_vol, calc_from_agen, not_supported };

//initialised in InitHists()
TH2D *fh_gen_nu;
TH2D *fh_gen_nub;
TH2D *fh_det_gandalf_nu;
TH2D *fh_det_shower_nu;
TH2D *fh_det_recolns_nu;
TH2D *fh_det_gandalf_nub;
TH2D *fh_det_shower_nub;
TH2D *fh_det_recolns_nub;

//gseagen and summary file parsers
GSGParser     *fG;
SummaryParser *fS;

//*********************************************************************

/**
 *  This routine fills effective mass histograms.
 *
 * \param  summary_file  Summary file with reconstructed neutrino events.
 * \param  gseagen_file  GSeaGen file where all of the generated MC events are.
 * \param  flavor       Neutrino flavor. 0 - electron, 1 - muon, 2 - tauon.
 * \param  int_type      Interaction type. 0 - NC, 1 - CC.
 * \param  en_low        MC neutrino energy start. Either 1 or 3.
 * \param  run_nr        MC run number.
 * \param  veff_option   Select effective volume calculation. 0 - interaction volume, 1 - detector can, 2 - from Agen
 *
 */
void EffectiveMass(TString summary_file, TString gseagen_file, 
		   Int_t flavor, Int_t int_type, Int_t en_low, Int_t run_nr,
		   Int_t veff_option = 1) {

  if (veff_option >= not_supported) {
    cout << "ERROR! EffectiveMass() veff_option " << veff_option 
	 << " not supported. Stopping." << endl;
    return;
  }

  TString out_name = ParseInputs(flavor, int_type, en_low, run_nr);
  if (out_name == "") return;

  InitHists();

  //------------------------------------------------------
  //load shared library, init datafile parsers
  //------------------------------------------------------
  gSystem->Load("../common_software/libcommonsoft.so");
  
  fS = new SummaryParser(summary_file);
  fG = new GSGParser(gseagen_file);

  //------------------------------------------------------
  //loops over summary events and fill 'detected' histograms
  //------------------------------------------------------
  FillDetected(veff_option);

  //------------------------------------------------------
  //calculate the quantity rho*Veff/Ngen, stored in h_gen_nu bins
  //different methods, depends on how the effective volume is calculated
  //------------------------------------------------------

  if (veff_option == calc_from_agen) { 
    FillGenerated_Agen(flavor, int_type);                
  }
  else {
    FillGenerated_MCevts(veff_option);
  }

  //------------------------------------------------------
  //calculate the effective mass by dividing with h_gen_nu, write to output
  //------------------------------------------------------

  TFile *fout = new TFile(out_name, "RECREATE");
  fh_gen_nu ->Write();
  fh_gen_nub->Write();
  fh_det_gandalf_nu->Write();
  fh_det_shower_nu ->Write();
  fh_det_recolns_nu->Write();
  fh_det_gandalf_nub->Write();
  fh_det_shower_nub ->Write();
  fh_det_recolns_nub->Write();
  fout->Close();

  //------------------------------------------------------
  //cleanup
  //------------------------------------------------------
  delete fh_gen_nu;
  delete fh_gen_nub;
  delete fh_det_gandalf_nu;
  delete fh_det_shower_nu;
  delete fh_det_recolns_nu;
  delete fh_det_gandalf_nub;
  delete fh_det_shower_nub;
  delete fh_det_recolns_nub;

  delete fS;
  delete fG;

}

//*********************************************************************

/** This routine parses the inputs to create the output file name. */
TString ParseInputs(Int_t flavor, Int_t int_type, Int_t en_low, Int_t run_nr) {

  map < Int_t, TString > f_to_name = { {0, "elec" }, 
				       {1, "muon" },
				       {2, "tau"  } };

  map < Int_t, TString > int_to_name = { {0, "NC" }, 
					 {1, "CC" } };

  map < Int_t, TString > en_to_name = { {1, "1-5GeV"   }, 
					{3, "3-100GeV" } };
  
  
  if ( f_to_name.find(flavor) == f_to_name.end() ) {
    cout << "ERROR! flavor " << flavor << " not supported." << endl;
    return "";
  }

  if ( int_to_name.find(int_type) == int_to_name.end() ) {
    cout << "ERROR! int_type " << int_type << " not supported." << endl;
    return "";
  }

  if ( en_to_name.find(en_low) == en_to_name.end() ) {
    cout << "ERROR! en_low " << en_low << " not supported." << endl;
    return "";
  }

  TString out_name = "output/effmass_" + f_to_name[flavor] + "-" + int_to_name[int_type] + "_" + 
    en_to_name[en_low] + "_" + (TString)to_string(run_nr) + ".root";

  return out_name;

}

//*********************************************************************

void InitHists() {

  fh_gen_nu  = new TH2D("Generated_nu" , "Generated_nu" , 100, 0, 100, 200, -1, 1);
  fh_gen_nub = (TH2D*)fh_gen_nu->Clone("Generated_nub");
  fh_gen_nub->SetNameTitle("Generated_nub","Generated_nub");

  fh_det_gandalf_nu  = (TH2D*)fh_gen_nu->Clone("Detected_gandalf_nu");
  fh_det_shower_nu   = (TH2D*)fh_gen_nu->Clone("Detected_shower_nu");
  fh_det_recolns_nu  = (TH2D*)fh_gen_nu->Clone("Detected_recolns_nu");
  fh_det_gandalf_nu->SetNameTitle("Detected_gandalf_nu","Detected_gandalf_nu");
  fh_det_shower_nu ->SetNameTitle("Detected_shower_nu" ,"Detected_shower_nu");
  fh_det_recolns_nu->SetNameTitle("Detected_recolns_nu","Detected_recolns_nu");

  fh_det_gandalf_nub = (TH2D*)fh_gen_nu->Clone("Detected_gandalf_nub");
  fh_det_shower_nub  = (TH2D*)fh_gen_nu->Clone("Detected_shower_nub");
  fh_det_recolns_nub = (TH2D*)fh_gen_nu->Clone("Detected_recolns_nub");
  fh_det_gandalf_nub->SetNameTitle("Detected_gandalf_nub","Detected_gandalf_nub");
  fh_det_shower_nub ->SetNameTitle("Detected_shower_nub" ,"Detected_shower_nub");
  fh_det_recolns_nub->SetNameTitle("Detected_recolns_nub","Detected_recolns_nub");

}

//*********************************************************************

Bool_t VertexInCan(Double_t vx, Double_t vy, Double_t vz, Double_t Rcan, Double_t zcan_min, Double_t zcan_max) {

  Double_t r_vtx = TMath::Sqrt(vx*vx + vy*vy);

  return ( ( r_vtx <= Rcan ) && ( vz >= zcan_min ) && ( vz <= zcan_max ) );
}

//*********************************************************************

void FillDetected(Int_t veff_option) {

  //loops over summary events and fill 'detected' histograms

  for (Int_t i = 0; i < fS->fChain->GetEntries(); i++) {

    fS->fChain->GetEntry(i);

    //if can volume is used as the effective volume skip all events with vertices outside the can
    if (veff_option == can_vol) {
      if ( !VertexInCan(fS->MC_pos_x, fS->MC_pos_y, fS->MC_pos_z, fG->Rcan, fG->Zcan_min, fG->Zcan_max) ) continue;
    }

    //gandalf
    if ((Bool_t)fS->gandalf_is_good) {

      if (fS->MC_type > 0){ fh_det_gandalf_nu ->Fill(fS->gandalf_energy_nu, fS->gandalf_dir_z ); }
      else                { fh_det_gandalf_nub->Fill(fS->gandalf_energy_nu, fS->gandalf_dir_z ); }

    }

    //shower
    if ((Bool_t)fS->shower_is_good ) {

      if (fS->MC_type > 0){ fh_det_shower_nu ->Fill(fS->shower_energy_nu , fS->shower_dir_z  ); }
      else                { fh_det_shower_nub->Fill(fS->shower_energy_nu , fS->shower_dir_z  ); }
      
    }

    //recolns
    if ((Bool_t)fS->recolns_is_good) {
      if (fS->MC_type > 0){ fh_det_recolns_nu ->Fill(fS->recolns_energy_nu, fS->recolns_dir_z ); }
      else                { fh_det_recolns_nub->Fill(fS->recolns_energy_nu, fS->recolns_dir_z ); }
    }

  }

}

//*********************************************************************

void FillGenerated_MCevts(Int_t veff_option) {

  //loop over generated MC events

  for (Int_t i = 0; i < fG->fChain->GetEntries(); i++) {
    fG->fChain->GetEntry(i);

    //if can volume is used as the effective volume skip all events with vertices outside the can
    if (veff_option == can_vol) {
      if ( !VertexInCan(fG->Neutrino_V1, fG->Neutrino_V2, fG->Neutrino_V3, fG->Rcan, fG->Zcan_min, fG->Zcan_max) ) continue;
    }

    if (fG->Neutrino_PdgCode > 0) { fh_gen_nu ->Fill(fG->Neutrino_E, fG->Neutrino_D3); }
    else                          { fh_gen_nub->Fill(fG->Neutrino_E, fG->Neutrino_D3); }
  }

  //multiply h_gen_nu by Veff * rho_seawater
  //after this step hDetected divided by h_gen_nu gives the effective mass hist

  Double_t scale = 0;
  Double_t rho   = fG->Rho_seawater;
  
  if (veff_option == interaction_vol) scale = fG->Vint * rho;
  else if (veff_option == can_vol) scale = fG->Vcan * rho;

  fh_gen_nu ->Scale( scale );
  fh_gen_nub->Scale( scale );

}

//*********************************************************************

//work in progress
void FillGenerated_Agen(Int_t flavor, Int_t int_type) {

  CScalculator cs;
  cs.SelectInteraction(flavor, (Bool_t)int_type, 0);

  //angular phase factor, assuming diffuse flux
  Double_t I_ct_norm = 2 * TMath::Pi() * (fG->Ct_max - fG->Ct_min);
  
  //energy phase factor
  Double_t I_E_norm = 0;

  if (fG->E_power >= 0) {
    cout << "ERROR! FillGenerated_Agen() unexpected E_power "
	 << fG->E_power << ", stopping." << endl;
    return;
  }
  else if (fG->E_power == -1) {
    I_E_norm = TMath::Log( fG->E_max/fG->E_min );
  }
  else {
    I_E_norm = ( TMath::Power(fG->E_max, fG->E_power + 1) - TMath::Power(fG->E_min, fG->E_power + 1) )/(fG->E_power + 1);
  }

  Double_t Ntot = fG->Ngen;

  for (Int_t Ebin = 1; Ebin <= fh_gen_nu->GetXaxis()->GetNbins(); Ebin++) {
    for (Int_t ctbin = 1; ctbin <= fh_gen_nu->GetYaxis()->GetNbins(); ctbin++) {

      Double_t emin  = fh_gen_nu->GetXaxis()->GetBinLowEdge(Ebin);
      Double_t emax  = fh_gen_nu->GetXaxis()->GetBinUpEdge(Ebin);
      Double_t ctmin = fh_gen_nu->GetYaxis()->GetBinLowEdge(ctbin);
      Double_t ctmax = fh_gen_nu->GetYaxis()->GetBinUpEdge(ctbin);

      //fraction of the angular phase space in the given bin
      Double_t ct_frac = 2 * TMath::Pi() * TMath::Abs(ctmax - ctmin)/I_ct_norm;

      //fraction of the energy phase space in the given bin
      Double_t e_frac = 0;

      if (fG->E_power == -1) {
	e_frac = TMath::Log( emax/emin )/I_E_norm;
      }
      else { 
	e_frac = ( TMath::Power(emax, fG->E_power + 1) - TMath::Power(emin, fG->E_power + 1) )/(fG->E_power + 1)/I_E_norm;
      }

      //calculate the fraction of the events in the specific bin, multiply by the cross-section
      // h_gen_nu ->SetBinContent(Ebin, ctbin, Ntot * ct_frac * e_frac);
      // h_gen_nub->SetBinContent(Ebin, ctbin, Ntot * ct_frac * e_frac);

    }
  }

}
