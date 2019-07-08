#include "EffMass.h"
#include "FileHeader.h"
#include "NMHUtils.h"

#include "TFile.h"
#include "TF1.h"
#include "TError.h"

#include<stdexcept>
#include<iostream>

using namespace std;

// define the dummyfile keyword
const TString EffMass::DUMMYFILE = "dummy";

/** Constructor for creating the file with effective mass histograms (use for writing).
    
    In this case all of the member histograms are initialised to NULL. The generated and selected histograms are to be set by the `EffMass::SetGenAndSelH` function, after which `EffMass::WriteToFile` is to be called to create the file. With this constructor the effective mass histograms are not initialised.

    \param emin    Minimum energy range (typically defined by the MC simulations, for ORCA usually 1 GeV)
    \param emax    Maximum energy range (typically defined by the MC simulations, for ORCA usually 100 GeV)
    \param ctmin   Minimum cos-theta range (usually -1)
    \param ctmax   Maximum cos-theta range (usually 1)
    \param bymin   Minimum bjroken-y range (usually 0)
    \param bymax   Maximum bjorken-y range (usually 1)
    \param datatag Tag to identify production, usually propaged from apps/data_sorting scripts.

*/
EffMass::EffMass(Double_t emin, Double_t emax, Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax, TString datatag) {

  fDataTag = datatag;
  fRhoSW = 0; // rho not used in this mode, set to 0

  fEmin  = emin;
  fEmax  = emax;
  fCtmin = ctmin;
  fCtmax = ctmax;
  fBymin = bymin;
  fBymax = bymax;

  // this constructor is meant for use in write-out mode; histgrams initalised to NULL and need to be set
  // by the user with SetGenAndSelH
  for (auto f: fFlavMap) {
    for (auto i: fIntMap) {
      for (auto p: fPolMap) {
	fhGen[f.first][i.first][p.first]  = NULL;
	fhSel[f.first][i.first][p.first]  = NULL;
	fhMeff[f.first][i.first][p.first] = NULL;
	fVgen[f.first][i.first][p.first]  = 0;
      }
    }
  }

}

//********************************************************************

/** Constructor for reading the file with effective mass histograms and providing effective mass functionality (use for reading).

    In this case the generated and selected histograms are read from the file and the effective mass historams are initiated, using the binning given through the constructor. Then the functions `EffMass::GetMeff...` can be used to get the effective mass value for a neutrino type at certain energy, cos-theta, bjorken-y.

    As a special case, the constructor can be initialised as `EffMass em(EffMass::DUMMYFILE,...)`, in which case a dummy value of \f$ \rho * 8e6 * ( 1 - e^{-0.5*energy} )\f$ is returned as the effective mass. This is useful for test applications.

    \param fname   File where the effective mass data have been written
    \param nebins  Number of energy bins to use
    \param nctbins Number of cos-theta bins to use
    \param nbybins Number of bjorken-y bins to use
    \param rho_sw  Seawater density

 */
EffMass::EffMass(TString fname, Int_t nebins, Int_t nctbins, Int_t nbybins, Double_t rho_sw) {

  fRhoSW = rho_sw;

  if (DUMMYFILE == fname) {
    CreateDummyData(nebins, nctbins, nbybins);
    cout << "NOTICE EffMass::EffMass() running in dummy mode, not suitable for physics analyses" << endl;
  }
  else ReadFromFile(fname);

  CreateMeffHists(nebins, nctbins, nbybins);

}

//********************************************************************

/** Destructor */
EffMass::~EffMass() {

  // delete dynamically allocated histograms
  for (auto f: fFlavMap) {
    for (auto i: fIntMap) {
      for (auto p: fPolMap) {

	if ( fhGen [f.first][i.first][p.first] ) delete fhGen [f.first][i.first][p.first] ;
	if ( fhSel [f.first][i.first][p.first] ) delete fhSel [f.first][i.first][p.first] ;
	if ( fhMeff[f.first][i.first][p.first] ) delete fhMeff[f.first][i.first][p.first] ;

	
      }
    }
  }

  // all stuff generated dynamically can be/has been added to fHeap, call delete on all of it
  TIter next(&fHeap);
  TObject *obj = NULL;
  while ( (obj = next() ) ) if (obj) delete obj;

}

//********************************************************************

/** Function to set the 'generated' and 'selected' histogram that are used for effective mass calculation.

    The effective mass is calculated as hsel/hgen * vgen * rho. Neither hsel nor hgen should be scaled, scaling is performed inside the class.

    \param flavor Neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param iscc   NC = 0, CC = 1
    \param isnb   nu = 0, nub = 1
    \param hgen   Histogram with 'generated' events, typically gSeaGen events inside the can
    \param hsel   Histogram with 'selected' events, typically events in PID output tree from ECAP
    \param vgen   Size of the generation volume in m3 (typically the can)
    
 */
void EffMass::SetGenAndSelH(Int_t flavor, Bool_t iscc, Bool_t isnb, TH3D* hgen, TH3D* hsel, Double_t vgen) {
  
  if (flavor > TAU) {
    throw std::invalid_argument("ERROR! EffMass::SetGenAndSelH() unknown flavor " + (string)to_string(flavor));
  }

  if ( hgen == NULL || hsel == NULL ) {
    throw std::invalid_argument("ERROR! EffMass::SetGenAndSelH() null pointer provided as argument.");
  }

  if ( hgen->GetEntries() < hsel->GetEntries() ) {
    throw std::invalid_argument("ERROR! EffMass::SetGenAndSelH() generated histogram has less entries than selected histogram.");
  }

  if ( !NMHUtils::BinsMatch(hgen, hsel) ) {
    throw std::invalid_argument("ERROR! EffMass::SetGenAndSelH() generated and selected histogram have different binning.");
  }

  if ( hsel->GetEntries() == 0 ) {
    cout << "WARNING! EffMass::SetGenAndSelH() selected histogram has 0 entries." << endl;
  }

  if ( hsel->GetEntries() == hgen->GetEntries() ) {
    cout << "WARNING! EffMass::SetGenAndSelH() selected histogram and generated histogram have the same number of entries." << endl;
  }

  if ( fhGen[flavor][(UInt_t)iscc][(UInt_t)isnb] != NULL || fhSel[flavor][(UInt_t)iscc][(UInt_t)isnb] != NULL ) {
    cout << "WARNING! EffMass::SetGenAndSelH() generated and selected histograms have already been set." << endl;
  }

  fhGen[flavor][(UInt_t)iscc][(UInt_t)isnb] = (TH3D*)hgen->Clone();
  fhSel[flavor][(UInt_t)iscc][(UInt_t)isnb] = (TH3D*)hsel->Clone();
  fhGen[flavor][(UInt_t)iscc][(UInt_t)isnb]->SetDirectory(0);
  fhSel[flavor][(UInt_t)iscc][(UInt_t)isnb]->SetDirectory(0);

  fVgen[flavor][(UInt_t)iscc][(UInt_t)isnb] = vgen;

}

//********************************************************************

/** Function to write the 'generated' and 'selected' histograms to a file.

    \param fname
 */
void EffMass::WriteToFile(TString fname) {

  FileHeader header("EffMass");
  header.AddParameter("emin" , (TString)to_string(fEmin)  );
  header.AddParameter("emax" , (TString)to_string(fEmax)  );
  header.AddParameter("ctmin", (TString)to_string(fCtmin) );
  header.AddParameter("ctmax", (TString)to_string(fCtmax) );
  header.AddParameter("bymin", (TString)to_string(fBymin) );
  header.AddParameter("bymax", (TString)to_string(fBymax) );
  header.AddParameter("datatag", fDataTag);
  
  TFile fout(fname, "RECREATE");

  // if any of the histograms un-initialized or no entries throw exception
  for (auto f: fFlavMap) {
    for (auto i: fIntMap) {
      for (auto p: fPolMap) {

	TString typestr = f.second + "_" + i.second + "_" + p.second;
	TString genname = "hgen_" + typestr;
	TString selname = "hsel_" + typestr;
	TString volname = "vol_" + typestr;

	header.AddParameter(volname, (TString)to_string(fVgen[f.first][i.first][p.first]) );

	TH3D* hgen = fhGen[f.first][i.first][p.first];
	TH3D* hsel = fhSel[f.first][i.first][p.first];

	if (hgen == NULL || hsel == NULL) {
	  throw std::invalid_argument("EffMass::WriteToFile() Un-initialised selected or generated histogram for " + (string)typestr + ". Use EffMass::SetGenAndSelH to set generated and selected histograms for each flavor, cc/nc, nu/nubar combination before calling EffMass::WriteToFile." );
	}

	hgen->SetNameTitle(genname, genname);
	hsel->SetNameTitle(selname, selname);

	hgen->Write();
	hsel->Write();

      }
    }
  }

  header.WriteHeader(&fout);
  fout.Close();

}

//********************************************************************

/** Function to read the 'generated' and 'selected' histograms from the input file.

    \param fname Input file name
 */
void EffMass::ReadFromFile(TString fname) {

  if ( !NMHUtils::FileExists(fname) ) {
    throw std::invalid_argument("EffMass::ReadFromFile() cannot find " + (string)fname);
  }

  // read header data

  FileHeader h1("EffMass");
  h1.ReadHeader(fname);
  
  fDataTag = h1.GetParameter("datatag");

  fEmin  = std::stod( (string)h1.GetParameter("emin") );
  fEmax  = std::stod( (string)h1.GetParameter("emax") );
  fCtmin = std::stod( (string)h1.GetParameter("ctmin") );
  fCtmax = std::stod( (string)h1.GetParameter("ctmax") );
  fBymin = std::stod( (string)h1.GetParameter("bymin") );
  fBymax = std::stod( (string)h1.GetParameter("bymax") );

  // read the generated and selected histograms from the input file  

  TFile fin(fname, "READ");

  for (auto f: fFlavMap) {
    for (auto i: fIntMap) {
      for (auto p: fPolMap) {

	TString typestr = f.second + "_" + i.second + "_" + p.second;
	TString genname = "hgen_" + typestr;
	TString selname = "hsel_" + typestr;
	TString volname = "vol_" + typestr;

	TH3D* hgen = (TH3D*)fin.Get(genname);
	TH3D* hsel = (TH3D*)fin.Get(selname);

	if (hgen == NULL || hsel == NULL) {
	  throw std::invalid_argument("EffMass::ReadFromFile() could not find selected or generated histogram for " + (string)typestr );
	}

	// read volume; clone histograms and detatch from input file
	fVgen[f.first][i.first][p.first] = std::stod( (string)h1.GetParameter(volname) );
	fhGen[f.first][i.first][p.first] = (TH3D*)hgen->Clone();
	fhSel[f.first][i.first][p.first] = (TH3D*)hsel->Clone();
	fhGen[f.first][i.first][p.first]->SetDirectory(0);
	fhSel[f.first][i.first][p.first]->SetDirectory(0);

      }
    }
  }

  fin.Close();

}

//********************************************************************

/** 
    This function sets up the `selected` and `generated` histograms for use in dummy mode.

    In this case the function `EffMass::GetMeff(...)` returns \f$ \rho * 8e6 * ( 1 - e^{-0.5*energy} )\f$ for the effective mass at a certain energy. This reaches a plateu of approx. 8e6 around 20 GeV and roughly mimics a turn-on region below that energy.

*/
void EffMass::CreateDummyData(Int_t nebins, Int_t nctbins, Int_t nbybins) {

  //---------------------------------------------------------------------
  // configure ranges, binning and volume for dummy data
  //---------------------------------------------------------------------

  fDataTag = "dummy";

  Double_t dummyvol = 8e6; // dummy volume is 8 MTon
  Double_t alpha = -0.5;   // parameter that defines the "turn-on" region.
  fEmin  =  1;             // dummy range is 1-100 GeV
  fEmax  =  1e2;
  fCtmin = -1;
  fCtmax =  1;
  fBymin =  0;
  fBymax =  1;

  TF1 func("func", "1 - TMath::Exp( [0]*x )", fEmin, fEmax);
  func.SetParameter(0, alpha);

  for (auto f: fFlavMap) {
    for (auto i: fIntMap) {
      for (auto p: fPolMap) {

	TString typestr = f.second + "_" + i.second + "_" + p.second;
	TString genname = "hgen_" + typestr;
	TString selname = "hsel_" + typestr;
	TString volname = "vol_" + typestr;

	//---------------------------------------------------------------------
	// create dummy histograms, such that gen/sel is always 1; then `EffMass::CreateMeffHists` 
	// will create effective mass histograms that return an effective mass of dummyvol * rho
	//---------------------------------------------------------------------

	fVgen[f.first][i.first][p.first] = dummyvol;
	fhGen[f.first][i.first][p.first] = new TH3D(genname, genname, nebins, fEmin, fEmax, 
						    nctbins, fCtmin, fCtmax, nbybins, fBymin, fBymax);
	fhSel[f.first][i.first][p.first] = new TH3D(selname, selname, nebins, fEmin, fEmax, 
						    nctbins, fCtmin, fCtmax, nbybins, fBymin, fBymax);
	fhGen[f.first][i.first][p.first]->SetDirectory(0);
	fhSel[f.first][i.first][p.first]->SetDirectory(0);

	for (Int_t xbin = 1; xbin <= fhGen[f.first][i.first][p.first]->GetXaxis()->GetNbins(); xbin++) {
	  for (Int_t ybin = 1; ybin <= fhGen[f.first][i.first][p.first]->GetYaxis()->GetNbins(); ybin++) {
	    for (Int_t zbin = 1; zbin <= fhGen[f.first][i.first][p.first]->GetZaxis()->GetNbins(); zbin++) {

	      Double_t E = fhGen[f.first][i.first][p.first]->GetXaxis()->GetBinCenter( xbin );
	      
	      fhSel[f.first][i.first][p.first]->SetBinContent(xbin, ybin, zbin, func.Eval(E) );
	      fhGen[f.first][i.first][p.first]->SetBinContent(xbin, ybin, zbin, 1.0);

	    }
	  }
	}

      }
    }
  }

}

//********************************************************************

/** Function to create the effective mass histograms from 'generated' and 'selected' histograms with desired binning.

    \param nebins   Number of energy bins
    \param nctbins  Number of cos-theta bins
    \param nbybins  Number of bjorken-y bins

 */
void EffMass::CreateMeffHists(Int_t nebins, Int_t nctbins, Int_t nbybins) {


  for (auto f: fFlavMap) {
    for (auto i: fIntMap) {
      for (auto p: fPolMap) {

	//---------------------------------------------------------------------
	// create the histogram name
	//---------------------------------------------------------------------
	TString typestr  = f.second + "_" + i.second + "_" + p.second;
	TString meffname = "hmeff_" + typestr;

	TH3D* hgen = fhGen[f.first][i.first][p.first];
	TH3D* hsel = fhSel[f.first][i.first][p.first];

	if ( !NMHUtils::BinsMatch(hgen, hsel) ) {
	  throw std::logic_error("ERROR! EffMass::CreateMeffHists() selected and generated histograms have different binning for type " + (string)typestr);
	}

	//---------------------------------------------------------------------
	// calculate rebinning
	//---------------------------------------------------------------------
	Int_t existing_ebins  = hgen->GetXaxis()->GetNbins();
	Int_t existing_ctbins = hgen->GetYaxis()->GetNbins();
	Int_t existing_bybins = hgen->GetZaxis()->GetNbins();
	
	Int_t rebinning_ebins  = existing_ebins /nebins;
	Int_t rebinning_ctbins = existing_ctbins/nctbins;
	Int_t rebinning_bybins = existing_bybins/nbybins;
	
	// check that rebinning is valid
	if ( existing_ebins % nebins != 0 ) {
	  throw std::invalid_argument( "ERROR! EffMass::CreateMeffHists() energy axis nbins=" + to_string(existing_ebins) + " cannot be rebinned to " + to_string( nebins ) + ". Change the requested binning or recreate the input file to EffMass with more bins." );
	}
	
	if ( existing_ctbins % nctbins != 0 ) {
	  throw std::invalid_argument( "ERROR! EffMass::CreateMeffHists() cos-theta axis nbins=" + to_string(existing_ctbins) + " cannot be rebinned to " + to_string( nctbins ) + ". Change the requested binning or recreate the input file to EffMass with more bins." );
	}

	if ( existing_bybins % nbybins != 0 ) {
	  throw std::invalid_argument( "ERROR! EffMass::CreateMeffHists() bjorken-y axis nbins=" + to_string(existing_bybins) + " cannot be rebinned to " + to_string( nbybins ) + ". Change the requested binning or recreate the input file to EffMass with more bins." );
	}
	  
	//---------------------------------------------------------------------
	// create the meff histogram by cloning the selected histogram, rebinning it and dividing it
	// by a clone of the generated histogram. Latter is cloned to keep the hgen with original
	// binning in memory
	//---------------------------------------------------------------------
	TH3D* meff = (TH3D*)hsel->Clone();
	meff->SetNameTitle(meffname, meffname);
	meff->SetDirectory(0);
	
	TH3D* genRB = (TH3D*)hgen->Clone();
	genRB->SetDirectory(0);

	meff ->Rebin3D(rebinning_ebins, rebinning_ctbins, rebinning_bybins);
	genRB->Rebin3D(rebinning_ebins, rebinning_ctbins, rebinning_bybins);
	meff->Divide(genRB);

	// apply scaling
	meff->Scale( fVgen[f.first][i.first][p.first] * fRhoSW );

	fhMeff[f.first][i.first][p.first] = meff;
	delete genRB;
	
      }
    }
  }
  
}

//********************************************************************

/** Get the effective mass for a neutrino type at given energy, cos-theta, bjorken-y.

    This function looks up the bin in the corresponding effective mass histogram and returns the effective mass in Tons. If interpolation is requested, a value corresponding the exact input coordinates (E,ct,by) is interpolated by using `TH3::Interpolate`. The coordinate that is outside the interpolation range (bin centers of the first and last bin on the corresponding axis), but still within the axis range (first bin low edge and last bin high edge ) is set to the corresponding bin center, such that interpolation can still be performed.

    In both cases (interpolated and not interpolated) an error is raised if the coordinate point is outside of the axis range.

    The performance in interpolated mode is expected to be (considerably) worse. Due to issues with `TH3->Interpolate`, it is recommended to have at least 2 bins on each axis.

    \param flavor       Neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param iscc         NC = 0, CC = 1
    \param isnb         nu = 0, nub = 1
    \param E_true       True neutrino energy
    \param Ct_true      True cos-theta
    \param By_true      True bjorken-y
    \param interpolate  Interpolate between bin centers using `TH3::Interpolate` to get the exact value.
    \return Effective mass in Tons
 */
Double_t EffMass::GetMeff(Int_t flavor, Bool_t iscc, Bool_t isnb, 
			  Double_t E_true, Double_t Ct_true, Double_t By_true, Bool_t interpolate) {

  if (flavor > TAU) {
    throw std::invalid_argument("ERROR! EffMass::GetMeff() unknown flavor " + to_string(flavor));
  }

  if ( !InCoveredRange(E_true, Ct_true, By_true) ) {
    throw std::invalid_argument("ERROR! EffMass::GetMeff() point (E, ct, by) = (" + to_string(E_true) + ", " + to_string(Ct_true) + ", " + to_string(By_true) + ") is outside the covered range" );
  }

  TH3D* meff = fhMeff[flavor][(UInt_t)iscc][(UInt_t)isnb];

  if (meff == NULL) {
    throw std::logic_error("ERROR! EffMass::GetMeff() effective mass histogram is not initialised.");
  }

  if ( interpolate ) {

    if ( meff->GetXaxis()->GetNbins() == 1 || meff->GetYaxis()->GetNbins() == 1 ) {
      throw std::invalid_argument("ERROR! EffMass::GetMeff() for interpolation at least 2 energy and 2 cos-theta bins are required.");
    }
    
    Double_t _e  = BringToInterpolRange( meff->GetXaxis(), E_true );
    Double_t _ct = BringToInterpolRange( meff->GetYaxis(), Ct_true );
    Double_t _by = BringToInterpolRange( meff->GetZaxis(), By_true );

    // ROOT inheritance scheme is im-perfect, `TH3->Interpolate` does not work if there is only one bjorken-y bin
    // as it does not manage to call TH2->Interpolate itself, hence I do this expensive operation here.
    if ( meff->GetZaxis()->GetNbins() == 1 ) {
      gErrorIgnoreLevel = kFatal; // because otherwise ROOT gives some warning about sumw2 that is irrelevant for interpolation
      TH2D* p = (TH2D*)meff->Project3D("yx");
      return p->Interpolate( _e, _ct );
      delete p;
      gErrorIgnoreLevel = kWarning;
    }
    else {
      return meff->Interpolate( _e, _ct, _by );
    }
    
  }
  else {
  
    Int_t ebin  = meff->GetXaxis()->FindBin(E_true);
    Int_t ctbin = meff->GetYaxis()->FindBin(Ct_true);
    Int_t bybin = meff->GetZaxis()->FindBin(By_true);

    return meff->GetBinContent(ebin, ctbin, bybin);

  }
  
}

//********************************************************************

/** Get the effective mass for a neutrino type in given energy, cos-theta, bjorken-y bin

    Returns the effective mass in Tons from the specified bin. This function avoids the bin lookup calls performed in `EffMass::GetMeff` and is thus slightly faster. It is meant for use where performance is important, but typically this means that the calling program needs some notion about which bin content to fetch.

    \param flavor       Neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param iscc         NC = 0, CC = 1
    \param isnb         nu = 0, nub = 1
    \param E_truebin    True neutrino energy bin
    \param Ct_truebin   True cos-theta bin
    \param By_truebin   True bjorken-y bin
    \return Effective mass in Tons
 */
Double_t EffMass::GetMeffBC(Int_t flavor, Bool_t iscc, Bool_t isnb, 
			    Int_t E_truebin, Int_t Ct_truebin, Int_t By_truebin) {

  if (flavor > TAU) {
    throw std::invalid_argument("ERROR! EffMass::GetMeffBC() unknown flavor " + to_string(flavor));
  }

  if ( fhMeff[flavor][(UInt_t)iscc][(UInt_t)isnb] == NULL ) {
    throw std::logic_error("ERROR! EffMass::GetMeffBC() effective mass histogram is not initialised.");
  }

  return fhMeff[flavor][(UInt_t)iscc][(UInt_t)isnb]->GetBinContent(E_truebin, Ct_truebin, By_truebin);

}

//********************************************************************

/** Private function to check that the requested energy, cos-theta and bjorken-y are in the range covered by the effective mass file
    \param E_true  energy
    \param Ct_true cos-theta
    \param By_true bjorken-y
    \return true if in covered range, false otherwise
 */
Bool_t EffMass::InCoveredRange(Double_t E_true, Double_t Ct_true, Double_t By_true) {

  return ( ( E_true  >= fEmin  && E_true  <= fEmax  ) && 
	   ( Ct_true >= fCtmin && Ct_true <= fCtmax ) && 
	   ( By_true >= fBymin && By_true <= fBymax ) );

}

//********************************************************************

/** Private function to bring the value to the interpolation range of `TH3::Interpolate`. 

    If the value is <= than the bin center of the first bin or >= than the bin center of the last bin, the return value is set to the corresponding bin center value
    \param  axis pointer to `TAxis` 
    \param  val  value
    \return      value such that the parameter is inside the interpolation range
 */
Double_t EffMass::BringToInterpolRange(TAxis *axis, Double_t val) {

  Double_t ret;
  
  Double_t min = axis->GetBinCenter(1);
  Double_t max = axis->GetBinCenter( axis->GetNbins() );

  if      (val <= min) { ret = min + 1e-10; }
  else if (val >= max) { ret = max - 1e-10; }
  else                 { ret = val;         }

  return ret;

}

//********************************************************************

/** Fetch the 3D histogram with effective mass data.

    \param flavor       Neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param iscc         NC = 0, CC = 1
    \param isnb         nu = 0, nub = 1
    \return TH3D histogram with effective mass data in Tons
 */
TH3D* EffMass::GetMeff3DH(Int_t flavor, Bool_t iscc, Bool_t isnb) {

  if (flavor > TAU) {
    throw std::invalid_argument("ERROR! EffMass::GetMeff3DH() unknown flavor " + to_string(flavor));
  }

  return fhMeff[flavor][(UInt_t)iscc][(UInt_t)isnb];

}

//********************************************************************

/** Function to plot the effective mass vs energy, averaged over up-going events.
   
    \param flavor       Neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param iscc         NC = 0, CC = 1
    \param isnb         nu = 0, nub = 1
    \return `TGraphErrors` object with the effective mass curve.
 */
TGraphErrors* EffMass::AverageUpgoing(Int_t flavor, Bool_t iscc, Bool_t isnb) {

  if (flavor > TAU) {
    throw std::invalid_argument("ERROR! EffMass::AverageUpgoing() unknown flavor " + to_string(flavor));
  }

  TString typestr = fFlavMap[flavor] + "_" + fIntMap[(Int_t)iscc] + "_" + fPolMap[(Int_t)isnb];

  // first, project generated and selected to 2D -> average over bjorken-y
  TH2D* hgen2D = (TH2D*)fhGen[flavor][(UInt_t)iscc][(UInt_t)isnb]->Project3D("yx");
  TH2D* hsel2D = (TH2D*)fhSel[flavor][(UInt_t)iscc][(UInt_t)isnb]->Project3D("yx");
  hgen2D->SetDirectory(0);  
  hsel2D->SetDirectory(0);  

  // then project to energy axis in the range -1 to 0 -> average over cos-theta
  Int_t ctbin_min = 1;
  Int_t ctbin_max = hgen2D->GetYaxis()->FindBin(0.) - 1;

  TH1D* hgen1D = hgen2D->ProjectionX(typestr + "_gen_px", ctbin_min, ctbin_max, "e");
  TH1D* hsel1D = hsel2D->ProjectionX(typestr + "_sel_px", ctbin_min, ctbin_max, "e");
  hgen1D->SetDirectory(0);  
  hsel1D->SetDirectory(0);  

  // rebin
  Int_t xrb = hgen1D->GetXaxis()->GetNbins()/fhMeff[flavor][(UInt_t)iscc][(UInt_t)isnb]->GetXaxis()->GetNbins();
  hgen1D->Rebin( xrb );
  hsel1D->Rebin( xrb );

  // do the divide --> convert to effective mass; put to graph
  hsel1D->Divide(hgen1D);
  hsel1D->Scale( fVgen[flavor][(UInt_t)iscc][(UInt_t)isnb] * fRhoSW );

  TGraphErrors *g = new TGraphErrors();
  g->SetNameTitle("meff1d_" + typestr, "meff1d_" + typestr);
  g->GetXaxis()->SetTitle("Energy [GeV]");
  g->GetYaxis()->SetTitle("Effective mass [Ton]");
  g->SetMarkerStyle(21);

  for (Int_t xb = 1; xb <= hsel1D->GetXaxis()->GetNbins(); xb++) {
  
    g->SetPoint(xb-1, hsel1D->GetXaxis()->GetBinCenter(xb), hsel1D->GetBinContent(xb) );
    g->SetPointError(xb-1, 0, hsel1D->GetBinError(xb) );

  }

  fHeap.Add(hgen2D);
  fHeap.Add(hsel2D);
  fHeap.Add(hgen1D);
  fHeap.Add(hsel1D);
  fHeap.Add(g);

  return g;

}

//********************************************************************

/** This function returns a pointer to TH1D with energy on the x-axis and eff mass on the y-axis.
    \param flavor       Neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param iscc         NC = 0, CC = 1
    \param isnb         nu = 0, nub = 1
    \param ct           cos-theta value at which the slice is drawn
    \param by           bjorken-y value at which the slice is drawn
    \return `TH1D` with the effective mass curve for given input
*/
TH1D* EffMass::GetSlice(Int_t flavor, Bool_t iscc, Bool_t isnb, Double_t ct, Double_t by) {

  if (flavor > TAU) {
    throw std::invalid_argument("ERROR! EffMass::GetSlice() unknown flavor " + to_string(flavor));
  }

  TH3D* hmeff = fhMeff[flavor][(UInt_t)iscc][(UInt_t)isnb];

  TString nametitle = "meff_ct=" + (TString)to_string(ct) + "_by=" + (TString)to_string(by);
  Int_t ybin = hmeff->GetYaxis()->FindBin( ct );
  Int_t zbin = hmeff->GetZaxis()->FindBin( by );

  TH1D* slice = hmeff->ProjectionX(nametitle, ybin, ybin, zbin, zbin);
  fHeap.Add( slice );

  return slice;

}

//********************************************************************

/** Get the histogram with 'generated' events for given neutrino type

    \param flavor       Neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param iscc         NC = 0, CC = 1
    \param isnb         nu = 0, nub = 1
    \return             Pointer to a TH3D histogram.
 */
TH3D* EffMass::GetGen(Int_t flavor, Bool_t iscc, Bool_t isnb) {

  if (flavor > TAU) {
    throw std::invalid_argument("ERROR! EffMass::GetGen() unknown flavor " + to_string(flavor));
  }

  return fhGen[flavor][(UInt_t)iscc][(UInt_t)isnb];

}

//********************************************************************

/** Get the histogram with 'selected' events for given neutrino type

    \param flavor       Neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param iscc         NC = 0, CC = 1
    \param isnb         nu = 0, nub = 1
    \return             Pointer to a TH3D histogram.
 */
TH3D* EffMass::GetSel(Int_t flavor, Bool_t iscc, Bool_t isnb) {

  if (flavor > TAU) {
    throw std::invalid_argument("ERROR! EffMass::GetSel() unknown flavor " + to_string(flavor));
  }

  return fhSel[flavor][(UInt_t)iscc][(UInt_t)isnb];

}
