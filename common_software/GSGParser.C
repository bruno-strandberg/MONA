#define GSGParser_cxx
#include "GSGParser.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include<sstream>
#include <sys/stat.h>

void GSGParser::Loop()
{
//   In a ROOT session, you can do:
//      root> .L GSGParser.C
//      root> GSGParser t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

//**************************************************************************

/**
 * This function initialises the gSeaGen file in .root format for reading.
 *
 * \param  fname   Input file name.
 * \return         True if init successful, false otherwise.
 */
Bool_t GSGParser::InitRootFile(TString fname) {

  Bool_t InitOK = true;

  TFile *f = new TFile(fname, "READ");
  TTree *t = NULL;
  TTree *h = NULL;
  
  if ( f->IsOpen() ) { 
    t = (TTree*)f->Get("Events");
    h = (TTree*)f->Get("Header");
  }
  else { 
    cout << "ERROR! GSGParser::InitRootFile() coult not open file " << fname << endl; 
    InitOK = false;
  }

  if (t != NULL && h != NULL) { 

    //init the Events tree for reading
    Init(t);
    
    //use the header reader class to get the number of events and the generation volume
    GSGHeaderReader hr(h);
    hr.GetEntry(0);

    fNgen = hr.NTot;

    //variables for interaction volume calculation
    fRho_seawater = hr.RhoSW;
    Double_t h_water      = hr.HSeaWater;
    Double_t rho_rock     = hr.RhoSR;
    Double_t h_rock       = hr.HRock;
    Double_t r_int        = hr.RVol;  //radius of interaction volume

    //interaction volume corresponding to pure water
    //Double_t rho_water    = 1;
    /* fVint = TMath::Pi() * r_int * r_int *  */
    /*    (h_water * Rho_seawater/rho_water + h_rock * rho_rock/rho_water); */
    
    //interaction volume corresponding to sea water
    fVint = TMath::Pi() * r_int * r_int * 
      (h_water + h_rock * rho_rock/fRho_seawater);

    //variables for can calculation
    fZcan_min = hr.Can1;
    fZcan_max = hr.Can2;
    fRcan     = hr.Can3;
    
    fVcan = TMath::Pi() * fRcan * fRcan * (fZcan_max  - fZcan_min);

    //variables for effective volume calculation based on generation area
    fAgen     = hr.Agen;
    fE_min    = hr.EvMin;
    fE_max    = hr.EvMax;
    fCt_min   = hr.CtMin;
    fCt_max   = hr.CtMax;
    fE_power  = -hr.Alpha;

  }
  else { 
    cout << "ERROR! GSGParser::InitRootFile() could not find tree(s) Events or Header in file " << fname << endl;
    InitOK = false;
  }

  return InitOK;

}

//**************************************************************************

/**
 * A function to loop over the gSeaGen file irrespective of the format.
 *
 * A for loop over the .evt (plain text) file is difficult, because it is not trivial to
 * get the number of events. It is easiest to while loop until no more events are found.
 * This is precisely what the function NextEvtEvent() does. The .root files can be looped
 * over by using the GetEntry(i) function or calling fChain->GetEntry(i) directly. However,
 * for using the parser generally it is much more convenient to have one function to call
 * to loop over the events irrespective of the format. This function serves that purpose.
 *
 * \return True if next event found, false otherwise.
 */
Bool_t GSGParser::NextEvent() {

  if (fIsRootFile) {
    fEntry++;
    return fChain->GetEntry(fEntry);
  }
  else {
    return NextEvtEvent();
  }

}

//**************************************************************************

/**
 * This function initialises the gSeaGen file in .evt format for reading.
 *
 * \param  fname   Input file name.
 * \return         True if init successful, false otherwise.
 */
#ifdef WAANET
Bool_t GSGParser::InitEvtFile(TString fname) {

  struct stat buf;
  if ( (stat(fname, &buf) != 0) ) {
    cout << "ERROR! GSGParser::InitEvtFile() coult not open file " << fname << endl;
    return false;
  }

  fEvtFile = new EventFile( (string)fname );

  //perform the same calculation as in InitRootFile()

  fNgen = stod( fEvtFile->header.get_field("genvol","numberOfEvents") );

  fRho_seawater     = 1.0397500; //hardcoded, value from a .root file header
  Double_t rho_rock = 2.65;      //hardcoded, value from a .root file header
  Double_t h_water = stod( fEvtFile->header.get_field("genvol","zmax") );
  Double_t h_rock  = stod( fEvtFile->header.get_field("genvol","zmin") );
  Double_t r_int   = stod( fEvtFile->header.get_field("genvol","r") );
    
  fVint = TMath::Pi() * r_int * r_int * 
    (h_water + h_rock * rho_rock/fRho_seawater);

  fZcan_min = stod( fEvtFile->header.get_field("can","zmin") );
  fZcan_max = stod( fEvtFile->header.get_field("can","zmax") );
  fRcan     = stod( fEvtFile->header.get_field("can","r") );
    
  fVcan = TMath::Pi() * fRcan * fRcan * (fZcan_max  - fZcan_min);

  fE_min    = stod( fEvtFile->header.get_field("cut_nu","Emin") );
  fE_max    = stod( fEvtFile->header.get_field("cut_nu","Emax") );
  fCt_min   = stod( fEvtFile->header.get_field("cut_nu","cosTmin") );
  fCt_max   = stod( fEvtFile->header.get_field("cut_nu","cosTmax") );
  fE_power  = stod( fEvtFile->header.get_field("spectrum","alpha") );

  //generation area is stored in w1, which is not loaded until first next()
  //instead, fetch it directly from the file
  ifstream in;
  string line;
  TString str;
  Double_t w1, w2, w3;

  in.open(fname);
  while ( getline (in, line) ) {
    if ( line.find("weights") != string::npos ) { 
      stringstream(line) >> str >> w1 >> w2 >> w3;
      break;
    }
  }
  in.close();

  fAgen = w1;

  return true;

}

//**************************************************************************

/**
 * This function loads the next event of the gSeaGen file in .evt format.
 *
 * Neutrino direction, vertex and energy are fetched from the .evt file in aanet format
 * and stored to the corresponding variables in the gSeaGen root format.
 *
 */
Bool_t GSGParser::NextEvtEvent() {
  
  Bool_t nextLoaded = fEvtFile->next();
  if (!nextLoaded) return nextLoaded;
  
  Neutrino_V1       = fEvtFile->evt.mc_trks[0].pos.x;
  Neutrino_V2       = fEvtFile->evt.mc_trks[0].pos.y;
  Neutrino_V3       = fEvtFile->evt.mc_trks[0].pos.z;
  Neutrino_D1       = fEvtFile->evt.mc_trks[0].dir.x;
  Neutrino_D2       = fEvtFile->evt.mc_trks[0].dir.y;
  Neutrino_D3       = fEvtFile->evt.mc_trks[0].dir.z;
  Neutrino_E        = fEvtFile->evt.mc_trks[0].E;
  Neutrino_PdgCode  = fEvtFile->evt.mc_trks[0].type;
  
  return nextLoaded;
}
#else
Bool_t GSGParser::InitEvtFile(TString fname) { return false; }
Bool_t GSGParser::NextEvtEvent()             { return false; }
#endif
