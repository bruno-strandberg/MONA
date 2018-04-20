#include "SummaryParser.h"
#include "GSGParser.h"
#include "TSystem.h"
#include "TH2.h"
#include "TVector3.h"
#include "TMath.h"
#include "NMHUtils.h"
#include <sys/stat.h>

//*********************************************************************
//functions

TString  ParseInputs(Int_t evt_type, Int_t int_type, Int_t en_low, Int_t run_nr);
void     InitHists();
Bool_t   VertexInVol(Double_t vx, Double_t vy, Double_t vz, Double_t Rcan, Double_t zcan_min, Double_t zcan_max);
void     FillDetected(Int_t veff_option, Double_t atmmu_cut);
void     FillGenerated(Int_t veff_option);

//*********************************************************************
//global histograms and variables

enum effvol_options { interaction_vol = 0, can_vol, custom_vol, not_supported };

//initialised in InitHists()
TH2D *fh_gen_nu;
TH2D *fh_gen_nub;
TH2D *fh_gen_scaled_nu;
TH2D *fh_gen_scaled_nub;

TH2D *fh_det_nu;
TH2D *fh_det_nub;
TH2D *fh_det_gandalf_nu;
TH2D *fh_det_shower_nu;
TH2D *fh_det_gandalf_nub;
TH2D *fh_det_shower_nub;

//gseagen and summary file parsers
GSGParser     *fG;
SummaryParser *fS;

//*********************************************************************

/**
 *  This routine fills effective mass histograms.
 *
 * \param  summary_file  Summary file with reconstructed neutrino events.
 * \param  gseagen_file  GSeaGen file where all of the generated MC events are.
 * \param  flavor        Neutrino flavor. 0 - electron, 1 - muon, 2 - tauon.
 * \param  int_type      Interaction type. 0 - NC, 1 - CC.
 * \param  en_low        MC neutrino energy start. Either 1 or 3.
 * \param  run_nr        MC run number.
 * \param  atmmu_cut     PID cut to reject atmospheric muons (0 - very strict, 1 - all events pass).
 * \param  veff_option   Select effective volume calculation. 0 - interaction volume, 1 - detector can, 2 - custom volume.
 * \param  rvol          Radius of custom volume.
 * \param  zmin_vol      Minimum z of custom volume.
 * \param  zmax_vol      Maximum z of custom volume.
 *
 */
void EffMhists(TString summary_file, TString gseagen_file, 
	       Int_t flavor, Int_t int_type, Int_t en_low, Int_t run_nr,
	       Double_t atmmu_cut = 1, Int_t veff_option = 1,
	       Double_t rvol = 0., Double_t zmin_vol = 0., Double_t zmax_vol = 0.) {

  if (veff_option >= not_supported) {
    cout << "ERROR! EffectiveMass() veff_option " << veff_option 
	 << " not supported. Stopping." << endl;
    return;
  }

  TString out_name = ParseInputs(flavor, int_type, en_low, run_nr);
  if (out_name == "") return;

  struct stat buf1, buf2;
  if ( (stat(summary_file, &buf1) != 0) || (stat(summary_file, &buf2) != 0) ) {
    cout << "ERROR! EffMhists() input file(s) missing." << endl;
    return;
  }

  InitHists();

  //------------------------------------------------------
  //load shared library, init datafile parsers
  //------------------------------------------------------
  gSystem->Load("$NMHDIR/common_software/libnmhsoft.so");
  
  fS = new SummaryParser(summary_file);
  fG = new GSGParser(gseagen_file);

  //------------------------------------------------------
  //if custom volume is used, overwrite the values stored in GSPParser
  //------------------------------------------------------

  if (veff_option == custom_vol) {
    fG->fRcan      = rvol;
    fG->fZcan_min  = zmin_vol;
    fG->fZcan_max  = zmax_vol;
    fG->fVcan      = TMath::Pi() * TMath::Power(fG->fRcan, 2) * (fG->fZcan_max - fG->fZcan_min);
  }
  
  //------------------------------------------------------
  //loops over summary events and fill 'detected' histograms
  //------------------------------------------------------
  FillDetected(veff_option, atmmu_cut);

  //------------------------------------------------------
  //fill the 'generated' histograms
  //------------------------------------------------------
  FillGenerated(veff_option);

  //------------------------------------------------------
  //write out the histograms. Division of det/gen has to be done 
  //later, this allows easy combining of the outputs.
  //------------------------------------------------------

  TFile *fout = new TFile(out_name, "RECREATE");
  fh_gen_nu ->Write();
  fh_gen_nub->Write();
  fh_gen_scaled_nu ->Write();
  fh_gen_scaled_nub->Write();
  fh_det_nu ->Write();
  fh_det_nub->Write();
  fh_det_gandalf_nu->Write();
  fh_det_shower_nu ->Write();
  fh_det_gandalf_nub->Write();
  fh_det_shower_nub ->Write();

  fout->Close();

  //------------------------------------------------------
  //cleanup
  //------------------------------------------------------
  delete fh_gen_nu;
  delete fh_gen_nub;
  delete fh_gen_scaled_nu;
  delete fh_gen_scaled_nub;
  delete fh_det_nu;
  delete fh_det_nub;
  delete fh_det_gandalf_nu;
  delete fh_det_shower_nu;
  delete fh_det_gandalf_nub;
  delete fh_det_shower_nub;
  
  delete fS;
  delete fG;

}

//*********************************************************************

// This routine parses the inputs to create the output file name.
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

  TString out_name = "output/EffMhists_" + f_to_name[flavor] + "-" + int_to_name[int_type] + "_" + 
    en_to_name[en_low] + "_" + (TString)to_string(run_nr) + ".root";

  return out_name;

}

//*********************************************************************

// This routine initialises the (global) histograms
void InitHists() {

  // 'generated' histograms
  vector<Double_t> low_edges = NMHUtils::GetLogBins(120, 1, 100);

  fh_gen_nu  = new TH2D("Generated_nu", "Generated_nu", 120, &low_edges[0], 200, -1, 1);
  fh_gen_nub        = (TH2D*)fh_gen_nu->Clone("Generated_nub");
  fh_gen_scaled_nu  = (TH2D*)fh_gen_nu->Clone("Generated_scaled_nu");
  fh_gen_scaled_nub = (TH2D*)fh_gen_nu->Clone("Generated_scaled_nub");
  fh_gen_nub       ->SetNameTitle("Generated_nub","Generated_nub");
  fh_gen_scaled_nu ->SetNameTitle("Generated_scaled_nu" ,"Generated_scaled_nu");
  fh_gen_scaled_nub->SetNameTitle("Generated_scaled_nub","Generated_scaled_nub");

  // 'detected' histograms

  fh_det_nu          = (TH2D*)fh_gen_nu->Clone("Detected_nu");
  fh_det_gandalf_nu  = (TH2D*)fh_gen_nu->Clone("Detected_gandalf_nu");
  fh_det_shower_nu   = (TH2D*)fh_gen_nu->Clone("Detected_shower_nu");
  fh_det_nu->SetNameTitle("Detected_nu" ,"Detected_nu");
  fh_det_gandalf_nu->SetNameTitle("Detected_gandalf_nu","Detected_gandalf_nu");
  fh_det_shower_nu ->SetNameTitle("Detected_shower_nu" ,"Detected_shower_nu");

  fh_det_nub         = (TH2D*)fh_gen_nu->Clone("Detected_nub");
  fh_det_gandalf_nub = (TH2D*)fh_gen_nu->Clone("Detected_gandalf_nub");
  fh_det_shower_nub  = (TH2D*)fh_gen_nu->Clone("Detected_shower_nub");
  fh_det_nub->SetNameTitle("Detected_nub","Detected_nub");
  fh_det_gandalf_nub->SetNameTitle("Detected_gandalf_nub","Detected_gandalf_nub");
  fh_det_shower_nub ->SetNameTitle("Detected_shower_nub" ,"Detected_shower_nub");
  
}

//*********************************************************************

// Function to check whether vertex is inside a volume
Bool_t VertexInVol(Double_t vx, Double_t vy, Double_t vz, Double_t R, Double_t z_min, Double_t z_max) {

  Double_t r_vtx = TMath::Sqrt(vx*vx + vy*vy);

  return ( ( r_vtx <= R ) && ( vz >= z_min ) && ( vz <= z_max ) );
}

//*********************************************************************

// Function that fills 'detected' histograms
void FillDetected(Int_t veff_option, Double_t atmmu_cut) {

  //loops over summary events and fill 'detected' histograms

  for (Int_t i = 0; i < fS->fChain->GetEntries(); i++) {

    fS->fChain->GetEntry(i);
    
    //if volume cut is used exclude events with vertices outside the volume cut
    if (veff_option == can_vol || veff_option == custom_vol) {
      if ( !VertexInVol(fS->MC_pos_x, fS->MC_pos_y, fS->MC_pos_z, fG->fRcan, fG->fZcan_min, fG->fZcan_max) ) continue;
    }

    //reject events that look like atmospheric muons
    if ( fS->PID_muon_score > atmmu_cut ) continue; 

    //if you wish to set PID cut, do so here
    
    //mc_truth
    if ( fS->MC_type > 0 ) { fh_det_nu ->Fill( fS->MC_energy, -fS->MC_dir_z ); }
    else                   { fh_det_nub->Fill( fS->MC_energy, -fS->MC_dir_z ); }

    //gandalf uses MC truth, but checks for reco 'quality cuts'
    if ((Bool_t)fS->gandalf_ql1) {

      if (fS->MC_type > 0) { fh_det_gandalf_nu ->Fill( fS->MC_energy, -fS->MC_dir_z ); }
      else                 { fh_det_gandalf_nub->Fill( fS->MC_energy, -fS->MC_dir_z ); }

    }

    //shower uses MC truth, but checks for reco 'quality cuts'
    if ((Bool_t)fS->shower_ql1 ) {

      if (fS->MC_type > 0) { fh_det_shower_nu ->Fill( fS->MC_energy, -fS->MC_dir_z ); }
      else                 { fh_det_shower_nub->Fill( fS->MC_energy, -fS->MC_dir_z ); }
      
    }
    
  }

}

//*********************************************************************

// Function that fills 'generated' histograms, generation vol either interaction volume or can.
void FillGenerated(Int_t veff_option) {

  //loop over generated MC events

  while ( fG->NextEvent() ) {
    
    //if volume cut is used exclude events with vertices outside the volume cut
    if (veff_option == can_vol || veff_option == custom_vol) {
      if ( !VertexInVol(fG->Neutrino_V1, fG->Neutrino_V2, fG->Neutrino_V3, fG->fRcan, fG->fZcan_min, fG->fZcan_max) ) continue;
    }

    if (fG->Neutrino_PdgCode > 0) { 
      fh_gen_nu ->Fill(fG->Neutrino_E, -fG->Neutrino_D3); 
      fh_gen_scaled_nu ->Fill(fG->Neutrino_E, -fG->Neutrino_D3); 
    }
    else { 
      fh_gen_nub->Fill(fG->Neutrino_E, -fG->Neutrino_D3);
      fh_gen_scaled_nub->Fill(fG->Neutrino_E, -fG->Neutrino_D3);
    }
  }

  //multiply h_gen_nu by Veff * rho_seawater
  //after this step hDetected divided by h_gen_nu gives the effective mass hist

  Double_t scale = 0;
  Double_t rho   = fG->fRho_seawater;
  
  if (veff_option == interaction_vol) scale = fG->fVint * rho;
  else if (veff_option == can_vol || veff_option == custom_vol) scale = fG->fVcan * rho;

  fh_gen_scaled_nu ->Scale( 1./scale );
  fh_gen_scaled_nub->Scale( 1./scale );

}
