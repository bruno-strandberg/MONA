//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 22 11:56:02 2018 by ROOT version 6.12/06
// from TTree Events/ 
// found on file: ../../tmp/gSeaGen_muon-CC_3-100GeV_18.root
//////////////////////////////////////////////////////////

#ifndef GSGParser_h
#define GSGParser_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "GSGHeaderReader.h"
#include "TMath.h"
#include <iostream>
#include <stdexcept>

//aanet class to read evt files, WAANET set in makefile
#ifdef WAANET
#include "EventFile.hh"
#define AANETEXISTS 1
#else
class EventFile;
#define AANETEXISTS 0
#endif

using namespace std;

// Header file for the classes stored in the TTree if any.
#include "TObject.h"

/** A ROOT-generated class to parse gSeaGen data in ROOT format.
 *
 *  This class is automatically generated by doing Events->MakeClass() after opening 
 *  a gSeaGen output in root format. There are some added variables for effective mass,
 *  can volume and effective area calculation. There are slight modifications to the constructor.
 *  The tree structure is described in the KM3NeT wiki under gSeaGen.
 *
 */

class GSGParser {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxTracks = 35;
   static constexpr Int_t kMaxEarthLeptons = 2;
   static constexpr Int_t kMaxSysWgt_WParam = 1;

   // Declaration of leaf types
 //GSeaEvent       *Events;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Int_t           iEvt;
   Double_t        PScale;
   Bool_t          VerInCan;
   string          TargetName;
   Double_t        TargetMass;
   Int_t           ScattId;
   Int_t           InterId;
   Double_t        Bx;
   Double_t        By;
   Int_t           NTracks;
   UInt_t          Neutrino_fUniqueID;
   UInt_t          Neutrino_fBits;
   UInt_t          Neutrino_Id;
   Double_t        Neutrino_V1;
   Double_t        Neutrino_V2;
   Double_t        Neutrino_V3;
   Double_t        Neutrino_D1;
   Double_t        Neutrino_D2;
   Double_t        Neutrino_D3;
   Double_t        Neutrino_E;
   Double_t        Neutrino_T;
   Int_t           Neutrino_PdgCode;
   Int_t           Neutrino_Type;
   UInt_t          PrimaryLepton_fUniqueID;
   UInt_t          PrimaryLepton_fBits;
   UInt_t          PrimaryLepton_Id;
   Double_t        PrimaryLepton_V1;
   Double_t        PrimaryLepton_V2;
   Double_t        PrimaryLepton_V3;
   Double_t        PrimaryLepton_D1;
   Double_t        PrimaryLepton_D2;
   Double_t        PrimaryLepton_D3;
   Double_t        PrimaryLepton_E;
   Double_t        PrimaryLepton_T;
   Int_t           PrimaryLepton_PdgCode;
   Int_t           PrimaryLepton_Type;
   Int_t           Tracks_;
   UInt_t          Tracks_fUniqueID[kMaxTracks];   //[Tracks_]
   UInt_t          Tracks_fBits[kMaxTracks];   //[Tracks_]
   UInt_t          Tracks_Id[kMaxTracks];   //[Tracks_]
   Double_t        Tracks_V1[kMaxTracks];   //[Tracks_]
   Double_t        Tracks_V2[kMaxTracks];   //[Tracks_]
   Double_t        Tracks_V3[kMaxTracks];   //[Tracks_]
   Double_t        Tracks_D1[kMaxTracks];   //[Tracks_]
   Double_t        Tracks_D2[kMaxTracks];   //[Tracks_]
   Double_t        Tracks_D3[kMaxTracks];   //[Tracks_]
   Double_t        Tracks_E[kMaxTracks];   //[Tracks_]
   Double_t        Tracks_T[kMaxTracks];   //[Tracks_]
   Int_t           Tracks_PdgCode[kMaxTracks];   //[Tracks_]
   Int_t           Tracks_Type[kMaxTracks];   //[Tracks_]
   Int_t           EarthLeptons_;
   UInt_t          EarthLeptons_fUniqueID[kMaxEarthLeptons];   //[EarthLeptons_]
   UInt_t          EarthLeptons_fBits[kMaxEarthLeptons];   //[EarthLeptons_]
   UInt_t          EarthLeptons_Id[kMaxEarthLeptons];   //[EarthLeptons_]
   Double_t        EarthLeptons_V1[kMaxEarthLeptons];   //[EarthLeptons_]
   Double_t        EarthLeptons_V2[kMaxEarthLeptons];   //[EarthLeptons_]
   Double_t        EarthLeptons_V3[kMaxEarthLeptons];   //[EarthLeptons_]
   Double_t        EarthLeptons_D1[kMaxEarthLeptons];   //[EarthLeptons_]
   Double_t        EarthLeptons_D2[kMaxEarthLeptons];   //[EarthLeptons_]
   Double_t        EarthLeptons_D3[kMaxEarthLeptons];   //[EarthLeptons_]
   Double_t        EarthLeptons_E[kMaxEarthLeptons];   //[EarthLeptons_]
   Double_t        EarthLeptons_T[kMaxEarthLeptons];   //[EarthLeptons_]
   Int_t           EarthLeptons_PdgCode[kMaxEarthLeptons];   //[EarthLeptons_]
   Int_t           EarthLeptons_Type[kMaxEarthLeptons];   //[EarthLeptons_]
   Double_t        GenWeight;
   Double_t        EvtWeight;
   UInt_t          SysWgt_fUniqueID;
   UInt_t          SysWgt_fBits;
   vector<double>  SysWgt_WList;
   Int_t           SysWgt_WParam_;
   UInt_t          SysWgt_WParam_fUniqueID[kMaxSysWgt_WParam];   //[SysWgt.WParam_]
   UInt_t          SysWgt_WParam_fBits[kMaxSysWgt_WParam];   //[SysWgt.WParam_]
   string          SysWgt_WParam_name[kMaxSysWgt_WParam];
   vector<double>  SysWgt_WParam_wght[kMaxSysWgt_WParam];
   Double_t        ColumnDepth;
   Double_t        PEarth;
   Double_t        XSecMean;
   Double_t        MJD;
   Double_t        LST;

   // List of branches
   TBranch        *b_Events_fUniqueID;   //!
   TBranch        *b_Events_fBits;   //!
   TBranch        *b_Events_iEvt;   //!
   TBranch        *b_Events_PScale;   //!
   TBranch        *b_Events_VerInCan;   //!
   TBranch        *b_Events_TargetName;   //!
   TBranch        *b_Events_TargetMass;   //!
   TBranch        *b_Events_ScattId;   //!
   TBranch        *b_Events_InterId;   //!
   TBranch        *b_Events_Bx;   //!
   TBranch        *b_Events_By;   //!
   TBranch        *b_Events_NTracks;   //!
   TBranch        *b_Events_Neutrino_fUniqueID;   //!
   TBranch        *b_Events_Neutrino_fBits;   //!
   TBranch        *b_Events_Neutrino_Id;   //!
   TBranch        *b_Events_Neutrino_V1;   //!
   TBranch        *b_Events_Neutrino_V2;   //!
   TBranch        *b_Events_Neutrino_V3;   //!
   TBranch        *b_Events_Neutrino_D1;   //!
   TBranch        *b_Events_Neutrino_D2;   //!
   TBranch        *b_Events_Neutrino_D3;   //!
   TBranch        *b_Events_Neutrino_E;   //!
   TBranch        *b_Events_Neutrino_T;   //!
   TBranch        *b_Events_Neutrino_PdgCode;   //!
   TBranch        *b_Events_Neutrino_Type;   //!
   TBranch        *b_Events_PrimaryLepton_fUniqueID;   //!
   TBranch        *b_Events_PrimaryLepton_fBits;   //!
   TBranch        *b_Events_PrimaryLepton_Id;   //!
   TBranch        *b_Events_PrimaryLepton_V1;   //!
   TBranch        *b_Events_PrimaryLepton_V2;   //!
   TBranch        *b_Events_PrimaryLepton_V3;   //!
   TBranch        *b_Events_PrimaryLepton_D1;   //!
   TBranch        *b_Events_PrimaryLepton_D2;   //!
   TBranch        *b_Events_PrimaryLepton_D3;   //!
   TBranch        *b_Events_PrimaryLepton_E;   //!
   TBranch        *b_Events_PrimaryLepton_T;   //!
   TBranch        *b_Events_PrimaryLepton_PdgCode;   //!
   TBranch        *b_Events_PrimaryLepton_Type;   //!
   TBranch        *b_Events_Tracks_;   //!
   TBranch        *b_Tracks_fUniqueID;   //!
   TBranch        *b_Tracks_fBits;   //!
   TBranch        *b_Tracks_Id;   //!
   TBranch        *b_Tracks_V1;   //!
   TBranch        *b_Tracks_V2;   //!
   TBranch        *b_Tracks_V3;   //!
   TBranch        *b_Tracks_D1;   //!
   TBranch        *b_Tracks_D2;   //!
   TBranch        *b_Tracks_D3;   //!
   TBranch        *b_Tracks_E;   //!
   TBranch        *b_Tracks_T;   //!
   TBranch        *b_Tracks_PdgCode;   //!
   TBranch        *b_Tracks_Type;   //!
   TBranch        *b_Events_EarthLeptons_;   //!
   TBranch        *b_EarthLeptons_fUniqueID;   //!
   TBranch        *b_EarthLeptons_fBits;   //!
   TBranch        *b_EarthLeptons_Id;   //!
   TBranch        *b_EarthLeptons_V1;   //!
   TBranch        *b_EarthLeptons_V2;   //!
   TBranch        *b_EarthLeptons_V3;   //!
   TBranch        *b_EarthLeptons_D1;   //!
   TBranch        *b_EarthLeptons_D2;   //!
   TBranch        *b_EarthLeptons_D3;   //!
   TBranch        *b_EarthLeptons_E;   //!
   TBranch        *b_EarthLeptons_T;   //!
   TBranch        *b_EarthLeptons_PdgCode;   //!
   TBranch        *b_EarthLeptons_Type;   //!
   TBranch        *b_Events_GenWeight;   //!
   TBranch        *b_Events_EvtWeight;   //!
   TBranch        *b_Events_SysWgt_fUniqueID;   //!
   TBranch        *b_Events_SysWgt_fBits;   //!
   TBranch        *b_Events_SysWgt_WList;   //!
   TBranch        *b_Events_SysWgt_WParam_;   //!
   TBranch        *b_SysWgt_WParam_fUniqueID;   //!
   TBranch        *b_SysWgt_WParam_fBits;   //!
   TBranch        *b_SysWgt_WParam_name;   //!
   TBranch        *b_SysWgt_WParam_wght;   //!
   TBranch        *b_Events_ColumnDepth;   //!
   TBranch        *b_Events_PEarth;   //!
   TBranch        *b_Events_XSecMean;   //!
   TBranch        *b_Events_MJD;   //!
   TBranch        *b_Events_LST;   //!

   GSGParser(TString fname);
   virtual ~GSGParser();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   //------------------------------------------------------
   //bstr: functions/members
   //------------------------------------------------------
   Double_t fNgen;         //!< number of generated events
   Double_t fVint;         //!< interaction volume
   Double_t fRho_seawater; //!< sea water density

   Double_t fVcan;         //!< can volume
   Double_t fRcan;         //!< can radius
   Double_t fZcan_min;     //!< can minimum
   Double_t fZcan_max;     //!< can maximum

   Double_t fAgen;         //!< generation area
   Double_t fE_min;        //!< minimum generated neutrion energy
   Double_t fE_max;        //!< maximum generated neutrion energy
   Double_t fCt_min;       //!< minimum generated neutrino direction
   Double_t fCt_max;       //!< maximum generated neutrino direction
   Double_t fE_power;      //!< power of the energy spectrum, E^{-E_power}

   EventFile *fEvtFile;    //<! pointer aanet .evt file reader
   Bool_t     fIsRootFile; //<! flag to indicate whether gsg file in .root or .evt format
   Long64_t   fEntry;      //<! global index used by NextEvent() function
   Bool_t     InitRootFile(TString fname);
   Bool_t     InitEvtFile(TString fname);
   Bool_t     NextEvtEvent();
   Bool_t     NextEvent();


};

#endif

#ifdef GSGParser_cxx

/**
 * Constructor.
 * \param  fname  An input gSeaGen file (either in root or evt format).
 */
GSGParser::GSGParser(TString fname) : fChain(0)
{

  Bool_t InitOK = false;

  fEvtFile    = NULL;
  fIsRootFile = true;
  fEntry      = -1;

  if      ( fname.Contains(".root") ) {
    InitOK = InitRootFile(fname);
    fIsRootFile = true;
  }
  else if ( fname.Contains(".evt")  ) {

    fIsRootFile = false;

    if (AANETEXISTS > 0) {
      InitOK = InitEvtFile(fname);
    }
    else {
      InitOK = false;
      cout << "ERROR! GSGParser::GSGParser() for .evt files compile against aanet." << endl;
    }
    
  }

  if (!InitOK) {
    throw std::invalid_argument( "ERROR! GSGParser::GSGParser() file " + fname + " not found or trying to parse .evt files without compilation against aanet." );
  }

}

GSGParser::~GSGParser()
{
   if (!fChain) return;
   if ( fChain->GetCurrentFile() ) fChain->GetCurrentFile()->Close();
}

Int_t GSGParser::GetEntry(Long64_t entry)
{

  if (fIsRootFile) {
    if (fChain) return fChain->GetEntry(entry);
    else return 0;
  }
  else {
    cout << "ERROR! GSGParser::GetEntry() can be used only with .root files, use LoadNextEvt() for parsing .evt files." << endl;
    return 0;
  }

}

Long64_t GSGParser::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void GSGParser::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_Events_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_Events_fBits);
   fChain->SetBranchAddress("iEvt", &iEvt, &b_Events_iEvt);
   fChain->SetBranchAddress("PScale", &PScale, &b_Events_PScale);
   fChain->SetBranchAddress("VerInCan", &VerInCan, &b_Events_VerInCan);
   fChain->SetBranchAddress("TargetName", &TargetName, &b_Events_TargetName);
   fChain->SetBranchAddress("TargetMass", &TargetMass, &b_Events_TargetMass);
   fChain->SetBranchAddress("ScattId", &ScattId, &b_Events_ScattId);
   fChain->SetBranchAddress("InterId", &InterId, &b_Events_InterId);
   fChain->SetBranchAddress("Bx", &Bx, &b_Events_Bx);
   fChain->SetBranchAddress("By", &By, &b_Events_By);
   fChain->SetBranchAddress("NTracks", &NTracks, &b_Events_NTracks);
   fChain->SetBranchAddress("Neutrino.fUniqueID", &Neutrino_fUniqueID, &b_Events_Neutrino_fUniqueID);
   fChain->SetBranchAddress("Neutrino.fBits", &Neutrino_fBits, &b_Events_Neutrino_fBits);
   fChain->SetBranchAddress("Neutrino.Id", &Neutrino_Id, &b_Events_Neutrino_Id);
   fChain->SetBranchAddress("Neutrino.V1", &Neutrino_V1, &b_Events_Neutrino_V1);
   fChain->SetBranchAddress("Neutrino.V2", &Neutrino_V2, &b_Events_Neutrino_V2);
   fChain->SetBranchAddress("Neutrino.V3", &Neutrino_V3, &b_Events_Neutrino_V3);
   fChain->SetBranchAddress("Neutrino.D1", &Neutrino_D1, &b_Events_Neutrino_D1);
   fChain->SetBranchAddress("Neutrino.D2", &Neutrino_D2, &b_Events_Neutrino_D2);
   fChain->SetBranchAddress("Neutrino.D3", &Neutrino_D3, &b_Events_Neutrino_D3);
   fChain->SetBranchAddress("Neutrino.E", &Neutrino_E, &b_Events_Neutrino_E);
   fChain->SetBranchAddress("Neutrino.T", &Neutrino_T, &b_Events_Neutrino_T);
   fChain->SetBranchAddress("Neutrino.PdgCode", &Neutrino_PdgCode, &b_Events_Neutrino_PdgCode);
   fChain->SetBranchAddress("Neutrino.Type", &Neutrino_Type, &b_Events_Neutrino_Type);
   fChain->SetBranchAddress("PrimaryLepton.fUniqueID", &PrimaryLepton_fUniqueID, &b_Events_PrimaryLepton_fUniqueID);
   fChain->SetBranchAddress("PrimaryLepton.fBits", &PrimaryLepton_fBits, &b_Events_PrimaryLepton_fBits);
   fChain->SetBranchAddress("PrimaryLepton.Id", &PrimaryLepton_Id, &b_Events_PrimaryLepton_Id);
   fChain->SetBranchAddress("PrimaryLepton.V1", &PrimaryLepton_V1, &b_Events_PrimaryLepton_V1);
   fChain->SetBranchAddress("PrimaryLepton.V2", &PrimaryLepton_V2, &b_Events_PrimaryLepton_V2);
   fChain->SetBranchAddress("PrimaryLepton.V3", &PrimaryLepton_V3, &b_Events_PrimaryLepton_V3);
   fChain->SetBranchAddress("PrimaryLepton.D1", &PrimaryLepton_D1, &b_Events_PrimaryLepton_D1);
   fChain->SetBranchAddress("PrimaryLepton.D2", &PrimaryLepton_D2, &b_Events_PrimaryLepton_D2);
   fChain->SetBranchAddress("PrimaryLepton.D3", &PrimaryLepton_D3, &b_Events_PrimaryLepton_D3);
   fChain->SetBranchAddress("PrimaryLepton.E", &PrimaryLepton_E, &b_Events_PrimaryLepton_E);
   fChain->SetBranchAddress("PrimaryLepton.T", &PrimaryLepton_T, &b_Events_PrimaryLepton_T);
   fChain->SetBranchAddress("PrimaryLepton.PdgCode", &PrimaryLepton_PdgCode, &b_Events_PrimaryLepton_PdgCode);
   fChain->SetBranchAddress("PrimaryLepton.Type", &PrimaryLepton_Type, &b_Events_PrimaryLepton_Type);
   fChain->SetBranchAddress("Tracks", &Tracks_, &b_Events_Tracks_);
   fChain->SetBranchAddress("Tracks.fUniqueID", Tracks_fUniqueID, &b_Tracks_fUniqueID);
   fChain->SetBranchAddress("Tracks.fBits", Tracks_fBits, &b_Tracks_fBits);
   fChain->SetBranchAddress("Tracks.Id", Tracks_Id, &b_Tracks_Id);
   fChain->SetBranchAddress("Tracks.V1", Tracks_V1, &b_Tracks_V1);
   fChain->SetBranchAddress("Tracks.V2", Tracks_V2, &b_Tracks_V2);
   fChain->SetBranchAddress("Tracks.V3", Tracks_V3, &b_Tracks_V3);
   fChain->SetBranchAddress("Tracks.D1", Tracks_D1, &b_Tracks_D1);
   fChain->SetBranchAddress("Tracks.D2", Tracks_D2, &b_Tracks_D2);
   fChain->SetBranchAddress("Tracks.D3", Tracks_D3, &b_Tracks_D3);
   fChain->SetBranchAddress("Tracks.E", Tracks_E, &b_Tracks_E);
   fChain->SetBranchAddress("Tracks.T", Tracks_T, &b_Tracks_T);
   fChain->SetBranchAddress("Tracks.PdgCode", Tracks_PdgCode, &b_Tracks_PdgCode);
   fChain->SetBranchAddress("Tracks.Type", Tracks_Type, &b_Tracks_Type);
   fChain->SetBranchAddress("EarthLeptons", &EarthLeptons_, &b_Events_EarthLeptons_);
   fChain->SetBranchAddress("EarthLeptons.fUniqueID", EarthLeptons_fUniqueID, &b_EarthLeptons_fUniqueID);
   fChain->SetBranchAddress("EarthLeptons.fBits", EarthLeptons_fBits, &b_EarthLeptons_fBits);
   fChain->SetBranchAddress("EarthLeptons.Id", EarthLeptons_Id, &b_EarthLeptons_Id);
   fChain->SetBranchAddress("EarthLeptons.V1", EarthLeptons_V1, &b_EarthLeptons_V1);
   fChain->SetBranchAddress("EarthLeptons.V2", EarthLeptons_V2, &b_EarthLeptons_V2);
   fChain->SetBranchAddress("EarthLeptons.V3", EarthLeptons_V3, &b_EarthLeptons_V3);
   fChain->SetBranchAddress("EarthLeptons.D1", EarthLeptons_D1, &b_EarthLeptons_D1);
   fChain->SetBranchAddress("EarthLeptons.D2", EarthLeptons_D2, &b_EarthLeptons_D2);
   fChain->SetBranchAddress("EarthLeptons.D3", EarthLeptons_D3, &b_EarthLeptons_D3);
   fChain->SetBranchAddress("EarthLeptons.E", EarthLeptons_E, &b_EarthLeptons_E);
   fChain->SetBranchAddress("EarthLeptons.T", EarthLeptons_T, &b_EarthLeptons_T);
   fChain->SetBranchAddress("EarthLeptons.PdgCode", EarthLeptons_PdgCode, &b_EarthLeptons_PdgCode);
   fChain->SetBranchAddress("EarthLeptons.Type", EarthLeptons_Type, &b_EarthLeptons_Type);
   fChain->SetBranchAddress("GenWeight", &GenWeight, &b_Events_GenWeight);
   fChain->SetBranchAddress("EvtWeight", &EvtWeight, &b_Events_EvtWeight);
   fChain->SetBranchAddress("SysWgt.fUniqueID", &SysWgt_fUniqueID, &b_Events_SysWgt_fUniqueID);
   fChain->SetBranchAddress("SysWgt.fBits", &SysWgt_fBits, &b_Events_SysWgt_fBits);
   fChain->SetBranchAddress("SysWgt.WList", &SysWgt_WList, &b_Events_SysWgt_WList);
   fChain->SetBranchAddress("SysWgt.WParam", &SysWgt_WParam_, &b_Events_SysWgt_WParam_);
   fChain->SetBranchAddress("SysWgt.WParam.fUniqueID", &SysWgt_WParam_fUniqueID, &b_SysWgt_WParam_fUniqueID);
   fChain->SetBranchAddress("SysWgt.WParam.fBits", &SysWgt_WParam_fBits, &b_SysWgt_WParam_fBits);
   fChain->SetBranchAddress("SysWgt.WParam.name", &SysWgt_WParam_name, &b_SysWgt_WParam_name);
   fChain->SetBranchAddress("SysWgt.WParam.wght", &SysWgt_WParam_wght, &b_SysWgt_WParam_wght);
   fChain->SetBranchAddress("ColumnDepth", &ColumnDepth, &b_Events_ColumnDepth);
   fChain->SetBranchAddress("PEarth", &PEarth, &b_Events_PEarth);
   fChain->SetBranchAddress("XSecMean", &XSecMean, &b_Events_XSecMean);
   fChain->SetBranchAddress("MJD", &MJD, &b_Events_MJD);
   fChain->SetBranchAddress("LST", &LST, &b_Events_LST);
   Notify();
}

Bool_t GSGParser::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GSGParser::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GSGParser::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GSGParser_cxx
