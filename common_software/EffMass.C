#include "EffMass.h"
#include "FileHeader.h"
#include "NMHUtils.h"

#include "TFile.h"

#include<stdexcept>
#include<iostream>

using namespace std;

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

    \param fname   File where the effective mass data have been written
    \param nebins  Number of energy bins to use
    \param nctbins Number of cos-theta bins to use
    \param nbybins Number of bjorken-y bins to use

 */
EffMass::EffMass(TString fname, Int_t nebins, Int_t nctbins, Int_t nbybins, Double_t rho_sw) {

  fRhoSW = rho_sw;
  ReadFromFile(fname);
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

    \fname Input file name
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

/** Function to create the effective mass histograms from 'generated' and 'selected' histograms with desired binning.

    \param nebins   Number of energy bins
    \param nctbins  Number of cos-theta bins
    \param nbybins  Number of bjorken-y bins
    \param rho_sw   Sea-water density

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

    This function looks up the bin in the corresponding effective mass histogram and returns the effective mass in Tons.

    \param flavor       Neutrino flavor (0 - elec, 1 - muon, 2 - tau)
    \param iscc         NC = 0, CC = 1
    \param isnb         nu = 0, nub = 1
    \param E_true       True neutrino energy
    \param Ct_true      True cos-theta
    \param By_true      True bjorken-y
    \param interpolate  Interpolate between bin centers to get the exact value (not supported yet)
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

  if (interpolate) {
    cout << "WARNING! EffMass::GetMeff() interpolation not yet implemented, returning content of corresponding bin." << endl;
  }

  TH3D* meff = fhMeff[flavor][(UInt_t)iscc][(UInt_t)isnb];

  if (meff == NULL) {
    throw std::logic_error("ERROR! EffMass::GetMeff() effective mass histogram is not initialised.");
  }
  
  Int_t ebin  = meff->GetXaxis()->FindBin(E_true);
  Int_t ctbin = meff->GetYaxis()->FindBin(Ct_true);
  Int_t bybin = meff->GetZaxis()->FindBin(By_true);

  return meff->GetBinContent(ebin, ctbin, bybin);

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

Bool_t EffMass::InCoveredRange(Double_t E_true, Double_t Ct_true, Double_t By_true) {

  return ( ( E_true  >= fEmin  && E_true  <= fEmax  ) && 
	   ( Ct_true >= fCtmin && Ct_true <= fCtmax ) && 
	   ( By_true >= fBymin && By_true <= fBymax ) );

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
    throw std::invalid_argument("ERROR! EffMass::AverageUpgoing() unknown flavor " + to_string(flavor));
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
