
// NMH headers
#include "SummaryParser.h"
#include "GSGParser.h"
#include "NMHUtils.h"
#include "FileHeader.h"

// root headers
#include "TH3.h"
#include "TVector3.h"
#include "TMath.h"

// jpp headers
#include "JTools/JRange.hh"
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

// cpp headers
#include <stdexcept>

using namespace JTOOLS;

namespace EFFMASS {

  //---------------------------------------------------------
  //functions
  //---------------------------------------------------------

  TString  CreateOutputName(TString output_dir, Int_t evt_type, Int_t int_type, Int_t en_low, Int_t en_high, Int_t run_nr);
  void     InitHists(Int_t ebins, Int_t ctbins, Int_t bybins, JRange<double> energy_range);
  Bool_t   VertexInVol(Double_t vx, Double_t vy, Double_t vz, Double_t Rcan, Double_t zcan_min, Double_t zcan_max);
  void     FillDetected(Int_t veff_option, Double_t atmmu_cut, Double_t noise_cut);
  void     FillGenerated(Int_t veff_option);

  //---------------------------------------------------------
  // histograms and variables
  //---------------------------------------------------------

  enum effvol_options { can_vol = 0, custom_vol, not_supported }; //!< effective volume options

  // cos-theta and bjorken-y always cover the full range; energy range comes from command line
  static const double fCtmin  =  -1.; //!< cos-theta minimum of histograms
  static const double fCtmax  =   1.; //!< cos-theta maximum of histograms
  static const double fBymin   =  0.; //!< bjorken-y minimum of histograms
  static const double fBymax   =  1.; //!< bjorken-y maximum of histograms

  // histograms, initialised in InitHists()
  TH3D *fh_gen_nu;           //!< histogram with gSeaGen events, nu
  TH3D *fh_gen_nub;          //!< histogram with gSeaGen events, nu-bar
  TH3D *fh_gen_scaled_nu;    //!< histogram with gSeaGen events, nu, scaled by 1/(V*rho)
  TH3D *fh_gen_scaled_nub;   //!< histogram with gSeaGen events, nub, scaled by 1/(V*rho)
  TH3D *fh_det_nu;           //!< histogram with selected (=summary) events, nu
  TH3D *fh_det_nub;          //!< histogram with selected (=summary) events, nub

  //gseagen and summary file parsers
  GSGParser     *fG;         //!< pointer to a gSeaGen file parser instance
  SummaryParser *fS;         //!< pointer to a `SummaryParser` instance

};

using namespace EFFMASS;

//*********************************************************************

/** This routine fills histograms for effective mass calculation.
    
    Effective mass is a trickly concept. It is defined as \f$ N_{\rm sel}/N_{\rm gen} \cdot V_{\rm gen} \rho_{\rm seawater} \f$. \f$ V_{\rm gen} \f$ is a generation volume in which the generated events are uniformly distributed. This is typically the can volume, granted that it is sufficiently large. \f$ N_{\rm gen} \f$ is then the number of generated events inside the generation volume and can be determined un-ambiguously from a `gSeaGen` file. However, anyone can define \f$ N_{\rm sel} \f$ with their favorite cuts.

    In this analysis, \f$ N_{\rm sel} \f$ is defined as the events that entered the summary PID tree distributed by ECAP. Noise cut, muon cut and further quality cuts are defined in `EventSelection` and `DetResponse` classes, which properly account for the reduction of events due to further cuts. <>Hence it is important to not cut on muons and noise, when generating effective mass histograms for use with, e.g., `FitUtil`!<\B>. The cuts are made possible only to allow to make plots compatible with, e.g., `Swim` and `ParamNMH`.

    The ambiguity would vanish, if we agreed collectively that we calculate the effective mass after the trigger. The events in PID summary tree have undergone some selection cuts that are poorly documented, but in principle their main function is to remove atmosperic muons with simple cuts and to check that reconstructions have worked in some basic way. If we calculated the effective mass after the trigger and had all of the data after the trigger in the summary file, all we would need to do in the NNMO package would be to use the new effective masses and take the selection cuts into account in `EventSelection`'s and `DetResponse`'s.
 
 */
int main(int argc, char **argv) {

  //------------------------------------------------------
  //parse command line arguments
  //------------------------------------------------------  
  TString        summary_file;
  TString        gseagen_file;
  TString        output_dir;
  JRange<double> energy_range;
  Int_t          nebins;
  Int_t          nctbins;
  Int_t          nbybins;
  Double_t       atmmu_cut;
  Double_t       noise_cut;
  Int_t          veff_option;
  Double_t       rvol;
  Double_t       zmin_vol;
  Double_t       zmax_vol;

  try {

    JParser<> zap("This routine fills histograms for effective mass calculation.");

    string rangecomment = "Energy range of the effective mass histogram. In ORCA simulations are usually done in two ranges, e.g. 1-5GeV and 3-100GeV. The range defined here needs to cover the full width of the two ranges, i.e. 1-100 GeV in this example.";
    string binningcomment = "This is chosen a large number that can be re-binned to a suitable number for analysis. The default value of 240 can be re-binned to 20, 24, 40, 48, ..., which are typically used in ORCA NMO.";

    zap['s'] = make_field(summary_file, "Summary file with reconstructed neutrino events");
    zap['g'] = make_field(gseagen_file, "GSeaGen file where all of the generated MC events are");
    zap['d'] = make_field(output_dir, "Directory where output file is written") = "output/";
    zap['E'] = make_field(nebins , "Number of energy bins. " + binningcomment) = 240;
    zap['C'] = make_field(nctbins, "Number of cos-theta bins. See comment for -E") = 240;
    zap['B'] = make_field(nbybins, "Number of bjorken-y bins. See comment for -E") = 24;
    zap['e'] = make_field(energy_range, rangecomment) = JRange<double>(1, 100);
    zap['a'] = make_field(atmmu_cut, "PID cut to reject atmospheric muons (0 - very strict, 1 - all events pass). NB! Be careful and read the documentation of the app!") = 1.;
    zap['n'] = make_field(noise_cut, "PID cut to reject noise-like events (0 - very strict, 1 - all events pass). NB! Be careful and read the documentation of the app!") = 1.;
    zap['v'] = make_field(veff_option, "Select effective volume calculation. 0 - detector can, 1 - custom volume.") = can_vol, custom_vol;
    zap['r'] = make_field(rvol , "Radius of custom volume") = 0.;
    zap['l'] = make_field(zmin_vol , "Minimum z of custom volume") = 0.;
    zap['u'] = make_field(zmax_vol , "Maximum z of custom volume") = 0.;

    zap(argc, argv);
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  if (veff_option >= not_supported) {
    throw std::invalid_argument("ERROR! EffMhists() veff_option " + to_string(veff_option) + " not supported." );
  }

  if ( !NMHUtils::FileExists(summary_file) || !NMHUtils::FileExists(gseagen_file) ) {
    throw std::invalid_argument("ERROR! EffMhists() input file(s) missing." );
  }

  //------------------------------------------------------
  //init datafile parsers, create output name and init hists
  //------------------------------------------------------  
  fS = new SummaryParser(summary_file);
  fG = new GSGParser(gseagen_file);

  // read the flavor and the interaction type from the first event in the summary file
  // this is not ideal, but for some reason the neutrino flavor and the interaction type (i.e. input settings)
  // are not stored in gSeaGen header.

  Int_t flavor   = fS->GetEvt(0)->Get_MC_type();
  Int_t int_type = fS->GetEvt(0)->Get_MC_is_CC();
  Int_t run_nr   = fG->fRunNr; // run number from gSeaGen
  Int_t en_low   = fG->fE_min; // minimum energy of the selected gSeaGen file file
  Int_t en_high  = fG->fE_max; // maximum energy of the selected gSeaGen file file

  if ( run_nr != fS->GetEvt(0)->Get_MC_runID() ) {
    throw std::invalid_argument("ERROR! EffMhists() gSeaGen and summary file run numbers are different.");
  }

  TString out_name = CreateOutputName(output_dir, flavor, int_type, en_low, en_high, run_nr);
  InitHists(nebins, nctbins, nbybins, energy_range);

  //------------------------------------------------------
  //if custom volume is used, overwrite the values stored in GSPParser (these are used in `FillGenerated`)
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
  FillDetected(veff_option, atmmu_cut, noise_cut);

  //------------------------------------------------------
  //fill the 'generated' histograms
  //------------------------------------------------------
  FillGenerated(veff_option);

  //------------------------------------------------------
  //write out the histograms. Division of det/gen has to be done 
  //later, this allows easy combining of the outputs.
  //------------------------------------------------------
  FileHeader h("EffMhists");
  h.ReadHeader(summary_file); // read the header from summary file, contains the tag
  h.AddParameter("summary_file", summary_file);
  h.AddParameter("gseagen_file", gseagen_file);
  h.AddParameter("emin", (TString)to_string(energy_range.getLowerLimit()) );
  h.AddParameter("emax", (TString)to_string(energy_range.getUpperLimit()) );
  h.AddParameter("atmmu_cut", (TString)to_string(atmmu_cut) );
  h.AddParameter("noise_cut", (TString)to_string(noise_cut) );
  h.AddParameter("veff_option", (TString)to_string(veff_option) );
  h.AddParameter("Rvol", (TString)to_string(fG->fRcan) );
  h.AddParameter("Zmin", (TString)to_string(fG->fZcan_min) );
  h.AddParameter("Zmax", (TString)to_string(fG->fZcan_max) );

  TFile *fout = new TFile(out_name, "RECREATE");
  fh_gen_nu ->Write();
  fh_gen_nub->Write();
  fh_gen_scaled_nu ->Write();
  fh_gen_scaled_nub->Write();
  fh_det_nu ->Write();
  fh_det_nub->Write();
  h.WriteHeader(fout);

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
  
  delete fS;
  delete fG;

}

//*********************************************************************

/**
   This function creates the output filename for given inputs.

   \param output_dir Directory where the output file is written
   \param nupdg      Neutrino pdg code
   \param int_type   Interaction type (0 = nc, 1 = cc)
   \param en_low     Start of energy range of gSeaGen simulation
   \param en_high    End of the energy range of gSeaGen simulation
   \param run_nr     gSeaGen run number

   \return          Output effective mass file name; returns empty if inputs not supported.

 */
TString EFFMASS::CreateOutputName(TString output_dir, Int_t nupdg, Int_t int_type, Int_t en_low, Int_t en_high, Int_t run_nr) {

  map < Int_t, TString > pdg_to_name = { {12, "elec" }, 
					 {14, "muon" },
					 {16, "tau"  } };

  map < Int_t, TString > int_to_name = { {0, "NC" }, 
					 {1, "CC" } };
  
  if ( pdg_to_name.find(TMath::Abs(nupdg)) == pdg_to_name.end() ) {
    throw std::invalid_argument("ERROR! EFFMASS::CreateOutputName() pdg code " + to_string(TMath::Abs(nupdg)) + " not supported.");
  }

  if ( int_to_name.find(int_type) == int_to_name.end() ) {
    throw std::invalid_argument("ERROR! EFFMASS::CreateOutputName() interaction type " + to_string(int_type) + " not supported.");
  }

  Int_t sysret = system("mkdir -p " + output_dir);

  if (sysret != 0) {
    throw std::logic_error("ERROR! EFFMASS::CreateOutputName() could not create dir " + (string)output_dir);
  }

  TString erange = to_string(en_low) + "-" + to_string(en_high) + "GeV";  

  TString out_name = output_dir + "/EffMhists_" + pdg_to_name[TMath::Abs(nupdg)] + "-" + int_to_name[int_type] + "_" + 
    erange + "_" + (TString)to_string(run_nr) + ".root";

  return out_name;

}

//*********************************************************************

/**
   Function that initialised the histograms used globally in this macro.
   \param ebins          Number of energy bins
   \param ctbins         Number of cos-theta bins
   \param bybins         Number of bjorken-y bins
   \param energy_range   Energy range of the Monte Carlo simulation
 */
void EFFMASS::InitHists(Int_t ebins, Int_t ctbins, Int_t bybins, JRange<double> energy_range) {

  // 'generated' histograms. Use fine binning here, as re-running this macro is time-consuming
  // suitable binning is chosen in the code that divides selected/generated

  Double_t emin   = energy_range.getLowerLimit();
  Double_t emax   = energy_range.getUpperLimit();
  
  vector<Double_t> e_edges  = NMHUtils::GetLogBins(ebins , emin  , emax);
  vector<Double_t> ct_edges = NMHUtils::GetBins   (ctbins, fCtmin, fCtmax);
  vector<Double_t> by_edges = NMHUtils::GetBins   (bybins, fBymin, fBymax);
  
  fh_gen_nu  = new TH3D("Generated_nu", "Generated_nu", 
			ebins, &e_edges[0], 
			ctbins, &ct_edges[0], 
			bybins, &by_edges[0]);

  fh_gen_nub        = (TH3D*)fh_gen_nu->Clone("Generated_nub");
  fh_gen_scaled_nu  = (TH3D*)fh_gen_nu->Clone("Generated_scaled_nu");
  fh_gen_scaled_nub = (TH3D*)fh_gen_nu->Clone("Generated_scaled_nub");
  fh_gen_nub       ->SetNameTitle("Generated_nub","Generated_nub");
  fh_gen_scaled_nu ->SetNameTitle("Generated_scaled_nu" ,"Generated_scaled_nu");
  fh_gen_scaled_nub->SetNameTitle("Generated_scaled_nub","Generated_scaled_nub");

  // 'detected' histograms
  fh_det_nu          = (TH3D*)fh_gen_nu->Clone("Detected_nu");
  fh_det_nu->SetNameTitle("Detected_nu" ,"Detected_nu");

  fh_det_nub         = (TH3D*)fh_gen_nu->Clone("Detected_nub");
  fh_det_nub->SetNameTitle("Detected_nub","Detected_nub");
  
}

//*********************************************************************

/**
   Function to check whether vertex is inside a volume

   \param vx     Vertex x
   \param vy     Vertex y
   \param vz     Vertex z
   \param R      Volume radius
   \param z_min  Minimum z coordinate of the volume
   \param z_max  Maximum z coordinate of the volume

   \return True if vertex inside volume, false otherwise.

 */
Bool_t EFFMASS::VertexInVol(Double_t vx, Double_t vy, Double_t vz, Double_t R, Double_t z_min, Double_t z_max) {

  Double_t r_vtx = TMath::Sqrt(vx*vx + vy*vy);

  return ( ( r_vtx <= R ) && ( vz >= z_min ) && ( vz <= z_max ) );
}

//*********************************************************************

/**
   Function that fills the 'detected' histograms of this macro.

   \param veff_option   Volume constraint option
   \param atmmu_cut     Atmospheric muon cut
   \param noise_cut     Noise cut
 */
void EFFMASS::FillDetected(Int_t veff_option, Double_t atmmu_cut, Double_t noise_cut) {

  //loops over summary events and fill 'detected' histograms

  for (Int_t i = 0; i < fS->GetTree()->GetEntries(); i++) {

    fS->GetTree()->GetEntry(i);
    SummaryEvent *evt = fS->GetEvt();

    // Do not cut on vertex of 'detected' events, this cannot be emulated in experimental data
    // still need to think about this...

    // if (veff_option == can_vol || veff_option == custom_vol) {
    //   TVector3 pos = evt->Get_MC_pos();
    //   if ( !VertexInVol(pos.x(), pos.y(), pos.z(), fG->fRcan, fG->fZcan_min, fG->fZcan_max) ) continue;
    // }

    //reject events that look like atmospheric muons (not used by default, accounted for in `EventSelection` and `DetResponse`)
    if ( evt->Get_RDF_muon_score() > atmmu_cut  ) continue; 
    if ( evt->Get_RDF_noise_score() > noise_cut ) continue; 

    //if you wish to set PID cut, do so here
    
    // fill nu and nubar 'detected' histograms
    if ( evt->Get_MC_type() > 0 ) { fh_det_nu ->Fill( evt->Get_MC_energy(), -evt->Get_MC_dir().z(), evt->Get_MC_bjorkeny() ); }
    else                          { fh_det_nub->Fill( evt->Get_MC_energy(), -evt->Get_MC_dir().z(), evt->Get_MC_bjorkeny() ); }
    
  }

}

//*********************************************************************

/**
   Function that fills the 'generated' histograms of this macro

   \param veff_option   Volume constraint option
 */
void EFFMASS::FillGenerated(Int_t veff_option) {

  //loop over generated MC events

  while ( fG->NextEvent() ) {
    
    //if volume cut is used exclude events with vertices outside the volume cut
    if (veff_option == can_vol || veff_option == custom_vol) {
      if ( !VertexInVol(fG->Neutrino_V1, fG->Neutrino_V2, fG->Neutrino_V3, fG->fRcan, fG->fZcan_min, fG->fZcan_max) ) continue;
    }

    if (fG->Neutrino_PdgCode > 0) { 
      fh_gen_nu ->Fill(fG->Neutrino_E, -fG->Neutrino_D3, fG->By ); 
      fh_gen_scaled_nu ->Fill(fG->Neutrino_E, -fG->Neutrino_D3, fG->By ); 
    }
    else { 
      fh_gen_nub->Fill(fG->Neutrino_E, -fG->Neutrino_D3, fG->By );
      fh_gen_scaled_nub->Fill(fG->Neutrino_E, -fG->Neutrino_D3, fG->By );
    }
  }

  //multiply h_gen_nu by Veff * rho_seawater
  //after this step hDetected divided by h_gen_nu gives the effective mass hist

  Double_t scale = 0;
  Double_t rho   = fG->fRho_seawater;
  
  if (veff_option == can_vol || veff_option == custom_vol) scale = fG->fVcan * rho;
  else { throw std::invalid_argument( "ERROR! FillGenerated() unknown volume option." ); }

  fh_gen_scaled_nu ->Scale( 1./scale );
  fh_gen_scaled_nub->Scale( 1./scale );

}
