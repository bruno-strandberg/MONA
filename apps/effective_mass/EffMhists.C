
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

#include "JGeometry3D/JTrack3E.hh"

#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

#include "evt/Evt.hh"
#include "evt/Trk.hh"
#include "evt/Head.hh"

#include "JAAnet/JHead.hh"
#include "JAAnet/JHeadToolkit.hh"
#include "JAAnet/JAAnetToolkit.hh"

#include "JSupport/JMultipleFileScanner.hh"
#include "JSupport/JMonteCarloFileSupportkit.hh"
#include "JSupport/JSupport.hh"

// cpp headers
#include <stdexcept>

using namespace JTOOLS;
using namespace JSUPPORT;
using namespace JPP;

namespace EFFMASS {

  //---------------------------------------------------------
  //functions
  //---------------------------------------------------------

  TString  CreateOutputName(TString output_dir, Int_t evt_type, Int_t int_type, Int_t en_low, Int_t en_high, Int_t run_nr);
  void     InitHists(Int_t ebins, Int_t ctbins, Int_t bybins, JRange<double> energy_range);
  Bool_t   VertexInVol(Double_t vx, Double_t vy, Double_t vz, Double_t Rcan, Double_t zcan_min, Double_t zcan_max);
  void     FillSelected(Int_t veff_option, Double_t atmmu_cut, Double_t noise_cut, Int_t flavor, Int_t int_type);
  void     FillGenerated(Int_t veff_option, Double_t rvol, Double_t zmin, Double_t zmax, Int_t flavor, Int_t int_type);

  //---------------------------------------------------------
  // histograms and variables
  //---------------------------------------------------------

  enum effvol_options { can_vol = 0, custom_vol, not_supported }; //!< effective volume options

  static const UInt_t GSG_CC = 2; //!< integer value to denote CC interaction in gSeaGen (see wiki)
  static const UInt_t GSG_NC = 3; //!< integer value to denote NC interaction in gSeaGen (see wiki)

  // cos-theta and bjorken-y always cover the full range; energy range comes from command line
  static const double fCtmin  =  -1.; //!< cos-theta minimum of histograms
  static const double fCtmax  =   1.; //!< cos-theta maximum of histograms
  static const double fBymin   =  0.; //!< bjorken-y minimum of histograms
  static const double fBymax   =  1.; //!< bjorken-y maximum of histograms

  // histograms, initialised in InitHists()
  TH3D *fh_gen_nu;           //!< histogram with gSeaGen events, nu
  TH3D *fh_gen_nub;          //!< histogram with gSeaGen events, nu-bar
  TH3D *fh_det_nu;           //!< histogram with selected (=summary) events, nu
  TH3D *fh_det_nub;          //!< histogram with selected (=summary) events, nub

  //gseagen and summary file parsers
  JMultipleFileScanner<Evt> fG; //!< gseagen scanner from JPP
  SummaryParser *fS;            //!< pointer to a `SummaryParser` instance

};

using namespace EFFMASS;

//*********************************************************************

/** This routine fills histograms for effective mass calculation.
    
    Effective mass is a trickly concept. It is defined as \f$ N_{\rm sel}/N_{\rm gen} \cdot V_{\rm gen} \rho_{\rm seawater} \f$. \f$ V_{\rm gen} \f$ is a generation volume in which the generated events are uniformly distributed. This is typically the can volume, granted that it is sufficiently large. \f$ N_{\rm gen} \f$ is then the number of generated events inside the generation volume and can be determined un-ambiguously from a `gSeaGen` file. However, anyone can define \f$ N_{\rm sel} \f$ with their favorite cuts.

    In this analysis, \f$ N_{\rm sel} \f$ is defined as the events that entered the summary PID tree distributed by ECAP. Noise cut, muon cut and further quality cuts are defined in `EventSelection` and `DetResponse` classes, which properly account for the reduction of events due to further cuts. <\B>Hence it is important to NOT cut on muons and noise, when generating effective mass histograms for use with, e.g., `FitUtil`!<\B>. The cuts are made possible only to allow to make plots compatible with, e.g., `Swim` and `ParamNMH`.

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
  fG = JMultipleFileScanner<Evt>((string)gseagen_file);
  const JHead aahead = getHeader(fG);

  // read the flavor and the interaction type from the first event in the summary file
  // this is not ideal, but for some reason the neutrino flavor and the interaction type (i.e. input settings)
  // are not stored in gSeaGen header.

  Int_t flavor      = fS->GetEvt(0)->Get_MC_type();
  Int_t int_type    = fS->GetEvt(0)->Get_MC_is_CC();
  Int_t run_nr      = aahead.start_run.run_id;        // run number from gSeaGen
  Double_t en_low   = aahead.cut_nu.Emin;             // minimum energy of the selected gSeaGen file file
  Double_t en_high  = aahead.cut_nu.Emax;             // maximum energy of the selected gSeaGen file file

  // if can volume is used over-write the vol dimensions with can dimensions
  if (veff_option == can_vol) {
    rvol     = aahead.can.r;
    zmin_vol = aahead.can.zmin;
    zmax_vol = aahead.can.zmax;
  }

  if ( run_nr != fS->GetEvt(0)->Get_MC_runID() ) {
    throw std::invalid_argument("ERROR! EffMhists() gSeaGen and summary file run numbers are different.");
  }

  if ( en_low != fS->GetEvt(0)->Get_MC_erange_start() ) {
    
    if ( flavor == 16 && TMath::Abs( en_low - fS->GetEvt(0)->Get_MC_erange_start() ) < 1 ) {
      cout << "WARNING! EffMhists() gSeaGen emin: " << en_low << ", summary emin: " 
	   << fS->GetEvt(0)->Get_MC_erange_start() << " for tau's; ignoring." << endl;
    }
    else {
      throw std::invalid_argument("ERROR! EffMhists() gSeaGen and summary file energy range mismatch.");
    }

  }

  TString out_name = CreateOutputName(output_dir, flavor, int_type, en_low, en_high, run_nr);
  InitHists(nebins, nctbins, nbybins, energy_range);
  
  //------------------------------------------------------
  //loops over summary events and fill 'selected' histograms
  //------------------------------------------------------
  FillSelected(veff_option, atmmu_cut, noise_cut, flavor, int_type);

  //------------------------------------------------------
  //fill the 'generated' histograms
  //------------------------------------------------------
  FillGenerated(veff_option, rvol, zmin_vol, zmax_vol, flavor, int_type);

  //------------------------------------------------------
  // write out the histograms. Division of det/gen and scaling has to be done 
  // later, this allows easy combining of the outputs
  //------------------------------------------------------
  FileHeader h("EffMhists");
  h.ReadHeader(summary_file); // read the header from summary file, contains the tag
  h.AddParameter("summary_file", summary_file);
  h.AddParameter("gseagen_file", gseagen_file);
  h.AddParameter("emin", (TString)to_string(energy_range.getLowerLimit()) );
  h.AddParameter("emax", (TString)to_string(energy_range.getUpperLimit()) );
  h.AddParameter("ctmin", (TString)to_string( fCtmin ) );
  h.AddParameter("ctmax", (TString)to_string( fCtmax ) );
  h.AddParameter("bymin", (TString)to_string( fBymin ) );
  h.AddParameter("bymax", (TString)to_string( fBymax ) );
  h.AddParameter("atmmu_cut", (TString)to_string(atmmu_cut) );
  h.AddParameter("noise_cut", (TString)to_string(noise_cut) );
  h.AddParameter("veff_option", (TString)to_string(veff_option) );
  h.AddParameter("rvol", (TString)to_string(rvol) );
  h.AddParameter("zmin", (TString)to_string(zmin_vol) );
  h.AddParameter("zmax", (TString)to_string(zmax_vol) );

  TFile *fout = new TFile(out_name, "RECREATE");
  fh_gen_nu ->Write();
  fh_gen_nub->Write();
  fh_det_nu ->Write();
  fh_det_nub->Write();
  h.WriteHeader(fout);

  fout->Close();

  //------------------------------------------------------
  //cleanup
  //------------------------------------------------------
  delete fh_gen_nu;
  delete fh_gen_nub;
  delete fh_det_nu;
  delete fh_det_nub;
  
  delete fS;

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
TString EFFMASS::CreateOutputName(TString output_dir, Int_t nupdg, Int_t int_type, 
				  Int_t en_low, Int_t en_high, Int_t run_nr) {

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

  TString out_name = output_dir + "/EffMhists_" + pdg_to_name[TMath::Abs(nupdg)] + "-" 
    + int_to_name[int_type] + "_" + erange + "_" + (TString)to_string(run_nr) + ".root";

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
  fh_gen_nub       ->SetNameTitle("Generated_nub","Generated_nub");

  // 'selected' histograms
  fh_det_nu          = (TH3D*)fh_gen_nu->Clone("Selected_nu");
  fh_det_nu->SetNameTitle("Selected_nu" ,"Selected_nu");

  fh_det_nub         = (TH3D*)fh_gen_nu->Clone("Selected_nub");
  fh_det_nub->SetNameTitle("Selected_nub","Selected_nub");
  
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
   Function that fills the 'selected' histograms of this macro.

   \param veff_option   Volume constraint option
   \param atmmu_cut     Atmospheric muon cut
   \param noise_cut     Noise cut
   \param flavor        nu flavor
   \param int_type      interaction type (cc/nc)
 */
void EFFMASS::FillSelected(Int_t veff_option, Double_t atmmu_cut, Double_t noise_cut, 
			   Int_t flavor, Int_t int_type) {

  //loops over summary events and fill 'selected' histograms

  for (Int_t i = 0; i < fS->GetTree()->GetEntries(); i++) {

    fS->GetTree()->GetEntry(i);
    SummaryEvent *evt = fS->GetEvt();

    // check that flavor and the interaction is always the same
    if ( TMath::Abs( flavor ) != TMath::Abs( evt->Get_MC_type() ) ) {
      throw std::invalid_argument("ERROR! FillSelected() flavor mismatch, expecting neutrinos of the same flavor (e.g. only 14, -14).");
    }

    if ( int_type != (Int_t)evt->Get_MC_is_CC() ) {
      throw std::invalid_argument("ERROR! FillSelected() interaction type mismatch, expecting only CC or NC events.");
    }

    // reject events that look like atmospheric muons; not used by default, accounted for in 
    // `EventSelection` and `DetResponse`, available to study their effects
    if ( evt->Get_RDF_muon_score() > atmmu_cut  ) continue; 
    if ( evt->Get_RDF_noise_score() > noise_cut ) continue; 
    
    // fill nu and nubar 'selected' histograms
    if ( evt->Get_MC_type() > 0 ) { 
      fh_det_nu ->Fill( evt->Get_MC_energy(), -evt->Get_MC_dir().z(), evt->Get_MC_bjorkeny() ); 
    }
    else { 
      fh_det_nub->Fill( evt->Get_MC_energy(), -evt->Get_MC_dir().z(), evt->Get_MC_bjorkeny() ); 
    }
    
  }

}

//*********************************************************************

/**
   Function that fills the 'generated' histograms of this macro

   \param veff_option   Volume constraint option
   \param rvol          Radius of the generation volume
   \param zmin          zmin of the generation volume
   \param zmax          zmax of the generation volume
   \param flavor        nu flavor
   \param int_type      interaction type (cc/nc)

 */
void EFFMASS::FillGenerated(Int_t veff_option, Double_t rvol, Double_t zmin, Double_t zmax,
			    Int_t flavor, Int_t int_type) {

  while ( fG.hasNext() ) {

    const Evt* event = fG.next();
    Trk nu_trk;
    
    if ( has_neutrino(*event) ) {
      nu_trk = get_neutrino(*event);
    }
    else {
      throw std::invalid_argument("ERROR! EFFMASS::FillGenerated() cannot find a neutrino in the MC event.");
    }

    // check that flavor and the interaction is always the same
    if ( TMath::Abs( flavor ) != TMath::Abs( nu_trk.type ) ) {
      throw std::invalid_argument("ERROR! FillGenerated() flavor mismatch, expecting neutrinos of same flavor (e.g. only 14, -14).");
    }

    Bool_t mismatch = ( ( ( int_type == 0 ) && ( (Int_t)nu_trk.getusr("cc") == GSG_CC ) ) ||
			( ( int_type == 1 ) && ( (Int_t)nu_trk.getusr("cc") == GSG_NC ) ) );

    if ( mismatch  ) {
      throw std::invalid_argument("ERROR! FillGenerated() interaction type mismatch, expecting only CC or NC events.");
    }

    //if volume cut is used exclude events with vertices outside the volume cut
    if (veff_option == can_vol || veff_option == custom_vol) {
      if ( !VertexInVol(nu_trk.pos.x, nu_trk.pos.y, nu_trk.pos.z, rvol, zmin, zmax) ) continue;
    }

    if (nu_trk.type > 0) { 
      fh_gen_nu        ->Fill( nu_trk.E, -nu_trk.dir.z, nu_trk.getusr("by") );
    }
    else { 
      fh_gen_nub       ->Fill( nu_trk.E, -nu_trk.dir.z, nu_trk.getusr("by") );
    }

  }

}
