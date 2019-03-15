
#include "EffMass.h"
#include "NMHUtils.h"
#include "FileHeader.h"

#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

#include "TMath.h"
#include "TFile.h"

#include <stdexcept>

/** Namespace that holds functions and variables for the `Combine` application. */
namespace COMBINE {

  /** Function to issue a system command and check for return.
      \param syscmd   System command
   */
  void SystemCmd( TString syscmd ) {
    Int_t sysret = system(syscmd);
    if (sysret != 0) {
      cout << "NOTICE COMBINE::SysretOK() returned " << sysret << " for cmd " << syscmd << endl;
    }
  }

  Double_t GetParameter(FileHeader &h, TString par_name);

  enum flavenum { ELEC = 0, MUON, TAU}; //!< flavor enumerator
  enum intenum  { NC = 0, CC};          //!< interaction types
  enum pols {NU=0, NUB};                //!< neutrino/anti-neutrino flags

};

/** This function operates on the outputs of `EffMhists` application. In practice, there are thousands of `gSeaGen` files and corresponding summary files, for each of those the `EffMhists` application is run, which creates an ouput (e.g. `EffMhists_elec-NC_3-100GeV_442.root`). This application takes the directory where the files are and merges the files by flavor. It will try to create files `Combined_elec-CC.root, Combined_muon-CC.root, Combined_tau-CC.root, Combined_elec-NC.root`. After this, it will try to use the `common_software/EffMass` class to create a single output file which can be used with other instances of the `common_software/EffMass` class to provide an effective mass calculator.

    For this application to work, some `EffMhists` outputs need to be available for all four neutrino interaction types (elec-CC, muon-CC, tau-CC, elec-NC). Otherwise the program will throw an exception.
*/
int main(const int argc, const char **argv) {

  using namespace COMBINE;

  //------------------------------------------------------
  //parse command line arguments
  //------------------------------------------------------  
  TString        input_dir;
  TString        output_dir;
  TString        output_file;
  Bool_t         use_datatag;

  try {

    JParser<> zap("This routine combines the outputs of EffMhists by flavor and interaction type and creates an output file for the class common_software/EffMass.");

    zap['i'] = make_field(input_dir      , "Directory with EffMhists outputs");
    zap['o'] = make_field(output_dir     , "Directory where combined outputs are written") = "combined_output/";
    zap['e'] = make_field(output_file    , "Name of the EffMass output file") = (TString)getenv("MONADIR") + "/data/eff_mass/EffMass.root";
    zap['t'] = make_field(use_datatag    , "Append the datatag to the output_file name. E.g. output/out.root becomes output/out_ORCA115_23x9m_ECAP0418.root (recommended)");

    if ( zap.read(argc, argv) != 0 ) return 1;
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  SystemCmd( "mkdir -p " + output_dir );

  //------------------------------------------------------
  // create combined output name and search string pairs for hadd
  //------------------------------------------------------  
  std::map<UInt_t, TString> flavors = { {ELEC, "elec"} , {MUON, "muon"} , {TAU, "tau"} };
  std::map<UInt_t, TString> ints    = { {NC, "NC"}, {CC, "CC"} };

  //tuples with combined output name, search string, flavor, interaction type
  enum members {COMBOUT=0, SEARCHSTR, FLAVOR, ISCC};
  vector< std::tuple<TString, TString, UInt_t, UInt_t> > P, Pe;

  for (auto f: flavors) {
    for (auto i: ints) {
      TString searchstr = input_dir + "/*" + f.second + "*" + i.second + "*";
      TString combout   = output_dir + "/Combined_" + f.second + "-" + i.second + ".root";
      P.push_back( std::make_tuple(combout, searchstr, f.first, i.first) );
    }
  }

  //------------------------------------------------------
  // call hadd
  //------------------------------------------------------  
  for (auto p: P) {

    TString syscmd;

    TString tmpfile = output_dir + "/tmp.txt";
    syscmd = "ls " + std::get<SEARCHSTR>(p) + " > " + tmpfile;
    SystemCmd ( syscmd );
    vector<TString> filelines = NMHUtils::ReadLines(tmpfile);
    SystemCmd ( "rm " + tmpfile );

    if (filelines.size() == 0) {
      cout << "NOTICE Combine() no files for selection " << std::get<SEARCHSTR>(p) << endl;
      continue;
    }
    else {
      cout << "NOTICE Combine() " << filelines.size() << " files for selection " << std::get<SEARCHSTR>(p) << endl;
    }

    syscmd = "hadd " + std::get<COMBOUT>(p) + " " + std::get<SEARCHSTR>(p);
    SystemCmd( syscmd );
    Pe.push_back(p); // vector with only existing files (this is because now we do not have muon, tau NC, but maybe in the future...)
  }

  if (Pe.size() == 0) {
    cout << "NOTICE Combine() no files created with hadd, returning." << endl;
    return 0;
  }

  //------------------------------------------------------
  // now use the effective mass class to create a combined output
  //------------------------------------------------------  

  // read in the header of the first file and get the limits
  FileHeader helper("helper");
  helper.ReadHeader( std::get<COMBOUT>(Pe[0]) );
  
  TString datatag = helper.GetParameter("datatag");

  Double_t emin  = GetParameter(helper, "emin");
  Double_t emax  = GetParameter(helper, "emax");
  Double_t ctmin = GetParameter(helper, "ctmin");
  Double_t ctmax = GetParameter(helper, "ctmax"); 
  Double_t bymin = GetParameter(helper, "bymin");
  Double_t bymax = GetParameter(helper, "bymax");

  Double_t atmmu_cut = GetParameter(helper, "atmmu_cut");
  Double_t noise_cut = GetParameter(helper, "noise_cut");

  EffMass em(emin, emax, ctmin, ctmax, bymin, bymax, datatag);

  //loop over the combined files and give histograms to `EffMass`

  for (auto p: Pe) {

    // check that MC ranges match
    FileHeader h("h");
    h.ReadHeader( std::get<COMBOUT>(p) );

    if ( emin  != GetParameter(h, "emin")  || emax  != GetParameter(h, "emax")  ||
	 ctmin != GetParameter(h, "ctmin") || ctmax != GetParameter(h, "ctmax") ||
	 bymin != GetParameter(h, "bymin") || bymax != GetParameter(h, "bymax") ) {
      throw std::logic_error("ERROR! Combine() MC range mismatch");
    }
    
    if ( atmmu_cut != GetParameter(h, "atmmu_cut") || noise_cut != GetParameter(h, "noise_cut") ) {
      throw std::logic_error("ERROR! Combine() muon or noise cut mismatch");
    }

    if ( datatag != h.GetParameter("datatag") ) {
      throw std::logic_error("ERROR! Combine() data tag mismatch; found " + (string)datatag + " and " + (string)h.GetParameter("datatag") );
    }

    //calculate the size of the generation volume (this can vary by flavor)
    Double_t rvol = GetParameter(h, "rvol");
    Double_t zmin = GetParameter(h, "zmin");
    Double_t zmax = GetParameter(h, "zmax");

    Double_t vgen = (zmax - zmin) * TMath::Pi() * rvol * rvol;

    TFile fin(std::get<COMBOUT>(p), "READ");

    if ( !fin.IsOpen() ) {
      cout << "NOTICE Combine() file " << fin.GetName() << " does not exist, continuing" << endl;
      continue;
    }

    TH3D *hgen_nu  = (TH3D*)fin.Get("Generated_nu");
    TH3D *hsel_nu  = (TH3D*)fin.Get("Selected_nu");
    TH3D *hgen_nub = (TH3D*)fin.Get("Generated_nub");
    TH3D *hsel_nub = (TH3D*)fin.Get("Selected_nub");

    em.SetGenAndSelH(std::get<FLAVOR>(p), std::get<ISCC>(p),  NU,  hgen_nu,  hsel_nu, vgen);
    em.SetGenAndSelH(std::get<FLAVOR>(p), std::get<ISCC>(p), NUB, hgen_nub, hsel_nub, vgen);
    
    // set elec-NC also for muon-NC and tau-NC
    if ( std::get<FLAVOR>(p) == ELEC && std::get<ISCC>(p) == NC ) {
      cout << "NOTICE Combine() setting elec NC eff masses also for muon and tau NC" << endl;

      em.SetGenAndSelH(MUON, NC,  NU,  hgen_nu,  hsel_nu, vgen);
      em.SetGenAndSelH(MUON, NC, NUB, hgen_nub, hsel_nub, vgen);
      em.SetGenAndSelH( TAU, NC,  NU,  hgen_nu,  hsel_nu, vgen);
      em.SetGenAndSelH( TAU, NC, NUB, hgen_nub, hsel_nub, vgen);
      
    }

    fin.Close();

  }

  //------------------------------------------------------
  // write the effective mass file
  //------------------------------------------------------  
  if (use_datatag) {
    output_file = output_file(0, output_file.Index(".root")) + "_" + datatag + ".root";
  }

  em.WriteToFile(output_file);

  // finally, add the muon and noise cuts to the output file for future reference
  FileHeader c("Combine");
  c.AddParameter(helper, "atmmu_cut");
  c.AddParameter(helper, "noise_cut");
  c.AddToFile(output_file);
}

//********************************************************************************************

/** This function looks up all values of the parameter in the header, checks that they are equal and returns the found value.
    \param h         Reference to a header where parameter is to be found
    \param par_name  Parameter name
*/
Double_t COMBINE::GetParameter(FileHeader &h, TString par_name) {

  auto parvalues = h.FindParValues(par_name);

  if (parvalues.size() == 0) {
    throw std::invalid_argument("ERROR! COMBINE::GetParameter() parameter " + (string)par_name + " cannot be found.");
  }

  Double_t par_value = std::stod( (string)parvalues[0].back() );

  for (auto p: parvalues) {
    Double_t pv = std::stod( (string)parvalues[0].back() );
    if ( pv != par_value ) {
      throw std::logic_error("ERROR! COMBINE::GetParameter() different values " + to_string(pv) + " and " + to_string(par_value) + " for parameter " + (string)par_name  );
    }
  }

  return par_value;

}
