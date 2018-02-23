//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 22 16:52:00 2018 by ROOT version 6.12/06
// from TTree Header/ 
// found on file: ../../tmp/gSeaGen_muon-CC_1-5GeV_10.root
//////////////////////////////////////////////////////////

#ifndef GSGHeaderReader_h
#define GSGHeaderReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TObject.h"
#include "vector"
#include "TNamed.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TAttMarker.h"
#include "TGraph.h"

#include<iostream>
using namespace std;

class GSGHeaderReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxSeaWaterComp = 10;
   static constexpr Int_t kMaxRockComp = 10;
   static constexpr Int_t kMaxMantleComp = 6;
   static constexpr Int_t kMaxCoreComp = 2;
   static constexpr Int_t kMaxFluxFiles = 2;

   // Declaration of leaf types
 //GenParam        *Header;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Int_t           RunNu;
   Int_t           RanSeed;
   string          InpXSecFile;
   Double_t        NTot;
   Double_t        EvMin;
   Double_t        EvMax;
   Double_t        CtMin;
   Double_t        CtMax;
   Double_t        Alpha;
   Int_t           NBin;
   Double_t        Can1;
   Double_t        Can2;
   Double_t        Can3;
   Double_t        HRock;
   Double_t        HSeaWater;
   Double_t        RVol;
   Double_t        NAbsLength;
   Double_t        AbsLength;
   Double_t        SiteDepth;
   Double_t        SiteLatitude;
   Double_t        SiteLongitude;
   Double_t        SeaBottomRadius;
   Double_t        GlobalGenWeight;
   Double_t        Agen;
   Double_t        RhoSW;
   Double_t        RhoSR;
   Double_t        TGen;
   string          PropCode;
   string          GenMode;
   string          SourceFile;
   Bool_t          GenMJD;
   string          SourceName;
   Double_t        Declination;
   Double_t        RightAscension;
   Double_t        SourceRadius;
   Double_t        MJDStart;
   Double_t        MJDStop;
   Int_t           SeaWaterComp_;
   Int_t           SeaWaterComp_first[kMaxSeaWaterComp];   //[SeaWaterComp_]
   Double_t        SeaWaterComp_second[kMaxSeaWaterComp];   //[SeaWaterComp_]
   Int_t           RockComp_;
   Int_t           RockComp_first[kMaxRockComp];   //[RockComp_]
   Double_t        RockComp_second[kMaxRockComp];   //[RockComp_]
   Int_t           MantleComp_;
   Int_t           MantleComp_first[kMaxMantleComp];   //[MantleComp_]
   Double_t        MantleComp_second[kMaxMantleComp];   //[MantleComp_]
   Int_t           CoreComp_;
   Int_t           CoreComp_first[kMaxCoreComp];   //[CoreComp_]
   Double_t        CoreComp_second[kMaxCoreComp];   //[CoreComp_]
   UInt_t          RangeSW_fUniqueID;
   UInt_t          RangeSW_fBits;
   TString         RangeSW_fName;
   TString         RangeSW_fTitle;
   Short_t         RangeSW_fLineColor;
   Short_t         RangeSW_fLineStyle;
   Short_t         RangeSW_fLineWidth;
   Short_t         RangeSW_fFillColor;
   Short_t         RangeSW_fFillStyle;
   Short_t         RangeSW_fMarkerColor;
   Short_t         RangeSW_fMarkerStyle;
   Float_t         RangeSW_fMarkerSize;
   Int_t           RangeSW_fNpoints;
   Double_t        RangeSW_fX[19];   //[RangeSW.fNpoints]
   Double_t        RangeSW_fY[19];   //[RangeSW.fNpoints]
   Double_t        RangeSW_fMinimum;
   Double_t        RangeSW_fMaximum;
   UInt_t          RangeSR_fUniqueID;
   UInt_t          RangeSR_fBits;
   TString         RangeSR_fName;
   TString         RangeSR_fTitle;
   Short_t         RangeSR_fLineColor;
   Short_t         RangeSR_fLineStyle;
   Short_t         RangeSR_fLineWidth;
   Short_t         RangeSR_fFillColor;
   Short_t         RangeSR_fFillStyle;
   Short_t         RangeSR_fMarkerColor;
   Short_t         RangeSR_fMarkerStyle;
   Float_t         RangeSR_fMarkerSize;
   Int_t           RangeSR_fNpoints;
   Double_t        RangeSR_fX[19];   //[RangeSR.fNpoints]
   Double_t        RangeSR_fY[19];   //[RangeSR.fNpoints]
   Double_t        RangeSR_fMinimum;
   Double_t        RangeSR_fMaximum;
   string          detector;
   vector<int>     NuList;
   Int_t           FluxFiles_;
   UInt_t          FluxFiles_fUniqueID[kMaxFluxFiles];   //[FluxFiles_]
   UInt_t          FluxFiles_fBits[kMaxFluxFiles];   //[FluxFiles_]
   Int_t           FluxFiles_NuType[kMaxFluxFiles];   //[FluxFiles_]
   string          FluxFiles_FluxSimul[kMaxFluxFiles];
   vector<string>  FluxFiles_FileNames[kMaxFluxFiles];
   Bool_t          PropMode;
   string          Drawing;
   string          gSeaGenVer;
   string          GenieVer;
   Long_t          RunTime;

   // List of branches
   TBranch        *b_Header_fUniqueID;   //!
   TBranch        *b_Header_fBits;   //!
   TBranch        *b_Header_RunNu;   //!
   TBranch        *b_Header_RanSeed;   //!
   TBranch        *b_Header_InpXSecFile;   //!
   TBranch        *b_Header_NTot;   //!
   TBranch        *b_Header_EvMin;   //!
   TBranch        *b_Header_EvMax;   //!
   TBranch        *b_Header_CtMin;   //!
   TBranch        *b_Header_CtMax;   //!
   TBranch        *b_Header_Alpha;   //!
   TBranch        *b_Header_NBin;   //!
   TBranch        *b_Header_Can1;   //!
   TBranch        *b_Header_Can2;   //!
   TBranch        *b_Header_Can3;   //!
   TBranch        *b_Header_HRock;   //!
   TBranch        *b_Header_HSeaWater;   //!
   TBranch        *b_Header_RVol;   //!
   TBranch        *b_Header_NAbsLength;   //!
   TBranch        *b_Header_AbsLength;   //!
   TBranch        *b_Header_SiteDepth;   //!
   TBranch        *b_Header_SiteLatitude;   //!
   TBranch        *b_Header_SiteLongitude;   //!
   TBranch        *b_Header_SeaBottomRadius;   //!
   TBranch        *b_Header_GlobalGenWeight;   //!
   TBranch        *b_Header_Agen;   //!
   TBranch        *b_Header_RhoSW;   //!
   TBranch        *b_Header_RhoSR;   //!
   TBranch        *b_Header_TGen;   //!
   TBranch        *b_Header_PropCode;   //!
   TBranch        *b_Header_GenMode;   //!
   TBranch        *b_Header_SourceFile;   //!
   TBranch        *b_Header_GenMJD;   //!
   TBranch        *b_Header_SourceName;   //!
   TBranch        *b_Header_Declination;   //!
   TBranch        *b_Header_RightAscension;   //!
   TBranch        *b_Header_SourceRadius;   //!
   TBranch        *b_Header_MJDStart;   //!
   TBranch        *b_Header_MJDStop;   //!
   TBranch        *b_Header_SeaWaterComp_;   //!
   TBranch        *b_SeaWaterComp_first;   //!
   TBranch        *b_SeaWaterComp_second;   //!
   TBranch        *b_Header_RockComp_;   //!
   TBranch        *b_RockComp_first;   //!
   TBranch        *b_RockComp_second;   //!
   TBranch        *b_Header_MantleComp_;   //!
   TBranch        *b_MantleComp_first;   //!
   TBranch        *b_MantleComp_second;   //!
   TBranch        *b_Header_CoreComp_;   //!
   TBranch        *b_CoreComp_first;   //!
   TBranch        *b_CoreComp_second;   //!
   TBranch        *b_Header_RangeSW_fUniqueID;   //!
   TBranch        *b_Header_RangeSW_fBits;   //!
   TBranch        *b_Header_RangeSW_fName;   //!
   TBranch        *b_Header_RangeSW_fTitle;   //!
   TBranch        *b_Header_RangeSW_fLineColor;   //!
   TBranch        *b_Header_RangeSW_fLineStyle;   //!
   TBranch        *b_Header_RangeSW_fLineWidth;   //!
   TBranch        *b_Header_RangeSW_fFillColor;   //!
   TBranch        *b_Header_RangeSW_fFillStyle;   //!
   TBranch        *b_Header_RangeSW_fMarkerColor;   //!
   TBranch        *b_Header_RangeSW_fMarkerStyle;   //!
   TBranch        *b_Header_RangeSW_fMarkerSize;   //!
   TBranch        *b_Header_RangeSW_fNpoints;   //!
   TBranch        *b_RangeSW_fX;   //!
   TBranch        *b_RangeSW_fY;   //!
   TBranch        *b_Header_RangeSW_fMinimum;   //!
   TBranch        *b_Header_RangeSW_fMaximum;   //!
   TBranch        *b_Header_RangeSR_fUniqueID;   //!
   TBranch        *b_Header_RangeSR_fBits;   //!
   TBranch        *b_Header_RangeSR_fName;   //!
   TBranch        *b_Header_RangeSR_fTitle;   //!
   TBranch        *b_Header_RangeSR_fLineColor;   //!
   TBranch        *b_Header_RangeSR_fLineStyle;   //!
   TBranch        *b_Header_RangeSR_fLineWidth;   //!
   TBranch        *b_Header_RangeSR_fFillColor;   //!
   TBranch        *b_Header_RangeSR_fFillStyle;   //!
   TBranch        *b_Header_RangeSR_fMarkerColor;   //!
   TBranch        *b_Header_RangeSR_fMarkerStyle;   //!
   TBranch        *b_Header_RangeSR_fMarkerSize;   //!
   TBranch        *b_Header_RangeSR_fNpoints;   //!
   TBranch        *b_RangeSR_fX;   //!
   TBranch        *b_RangeSR_fY;   //!
   TBranch        *b_Header_RangeSR_fMinimum;   //!
   TBranch        *b_Header_RangeSR_fMaximum;   //!
   TBranch        *b_Header_detector;   //!
   TBranch        *b_Header_NuList;   //!
   TBranch        *b_Header_FluxFiles_;   //!
   TBranch        *b_FluxFiles_fUniqueID;   //!
   TBranch        *b_FluxFiles_fBits;   //!
   TBranch        *b_FluxFiles_NuType;   //!
   TBranch        *b_FluxFiles_FluxSimul;   //!
   TBranch        *b_FluxFiles_FileNames;   //!
   TBranch        *b_Header_PropMode;   //!
   TBranch        *b_Header_Drawing;   //!
   TBranch        *b_Header_gSeaGenVer;   //!
   TBranch        *b_Header_GenieVer;   //!
   TBranch        *b_Header_RunTime;   //!

   GSGHeaderReader(TTree *tree=0);
   virtual ~GSGHeaderReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef GSGHeaderReader_cxx
GSGHeaderReader::GSGHeaderReader(TTree *tree) : fChain(0) 
{

   if (tree == 0) {
     cout << "ERROR! GSGHeaderReader::GSGHeaderReader() null pointer to the tree." << endl;
   }
   else {
     Init(tree);
   }
}

GSGHeaderReader::~GSGHeaderReader()
{
   if (!fChain) return;
   //delete fChain->GetCurrentFile();
}

Int_t GSGHeaderReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GSGHeaderReader::LoadTree(Long64_t entry)
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

void GSGHeaderReader::Init(TTree *tree)
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

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_Header_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_Header_fBits);
   fChain->SetBranchAddress("RunNu", &RunNu, &b_Header_RunNu);
   fChain->SetBranchAddress("RanSeed", &RanSeed, &b_Header_RanSeed);
   fChain->SetBranchAddress("InpXSecFile", &InpXSecFile, &b_Header_InpXSecFile);
   fChain->SetBranchAddress("NTot", &NTot, &b_Header_NTot);
   fChain->SetBranchAddress("EvMin", &EvMin, &b_Header_EvMin);
   fChain->SetBranchAddress("EvMax", &EvMax, &b_Header_EvMax);
   fChain->SetBranchAddress("CtMin", &CtMin, &b_Header_CtMin);
   fChain->SetBranchAddress("CtMax", &CtMax, &b_Header_CtMax);
   fChain->SetBranchAddress("Alpha", &Alpha, &b_Header_Alpha);
   fChain->SetBranchAddress("NBin", &NBin, &b_Header_NBin);
   fChain->SetBranchAddress("Can1", &Can1, &b_Header_Can1);
   fChain->SetBranchAddress("Can2", &Can2, &b_Header_Can2);
   fChain->SetBranchAddress("Can3", &Can3, &b_Header_Can3);
   fChain->SetBranchAddress("HRock", &HRock, &b_Header_HRock);
   fChain->SetBranchAddress("HSeaWater", &HSeaWater, &b_Header_HSeaWater);
   fChain->SetBranchAddress("RVol", &RVol, &b_Header_RVol);
   fChain->SetBranchAddress("NAbsLength", &NAbsLength, &b_Header_NAbsLength);
   fChain->SetBranchAddress("AbsLength", &AbsLength, &b_Header_AbsLength);
   fChain->SetBranchAddress("SiteDepth", &SiteDepth, &b_Header_SiteDepth);
   fChain->SetBranchAddress("SiteLatitude", &SiteLatitude, &b_Header_SiteLatitude);
   fChain->SetBranchAddress("SiteLongitude", &SiteLongitude, &b_Header_SiteLongitude);
   fChain->SetBranchAddress("SeaBottomRadius", &SeaBottomRadius, &b_Header_SeaBottomRadius);
   fChain->SetBranchAddress("GlobalGenWeight", &GlobalGenWeight, &b_Header_GlobalGenWeight);
   fChain->SetBranchAddress("Agen", &Agen, &b_Header_Agen);
   fChain->SetBranchAddress("RhoSW", &RhoSW, &b_Header_RhoSW);
   fChain->SetBranchAddress("RhoSR", &RhoSR, &b_Header_RhoSR);
   fChain->SetBranchAddress("TGen", &TGen, &b_Header_TGen);
   fChain->SetBranchAddress("PropCode", &PropCode, &b_Header_PropCode);
   fChain->SetBranchAddress("GenMode", &GenMode, &b_Header_GenMode);
   fChain->SetBranchAddress("SourceFile", &SourceFile, &b_Header_SourceFile);
   fChain->SetBranchAddress("GenMJD", &GenMJD, &b_Header_GenMJD);
   fChain->SetBranchAddress("SourceName", &SourceName, &b_Header_SourceName);
   fChain->SetBranchAddress("Declination", &Declination, &b_Header_Declination);
   fChain->SetBranchAddress("RightAscension", &RightAscension, &b_Header_RightAscension);
   fChain->SetBranchAddress("SourceRadius", &SourceRadius, &b_Header_SourceRadius);
   fChain->SetBranchAddress("MJDStart", &MJDStart, &b_Header_MJDStart);
   fChain->SetBranchAddress("MJDStop", &MJDStop, &b_Header_MJDStop);
   fChain->SetBranchAddress("SeaWaterComp", &SeaWaterComp_, &b_Header_SeaWaterComp_);
   fChain->SetBranchAddress("SeaWaterComp.first", SeaWaterComp_first, &b_SeaWaterComp_first);
   fChain->SetBranchAddress("SeaWaterComp.second", SeaWaterComp_second, &b_SeaWaterComp_second);
   fChain->SetBranchAddress("RockComp", &RockComp_, &b_Header_RockComp_);
   fChain->SetBranchAddress("RockComp.first", RockComp_first, &b_RockComp_first);
   fChain->SetBranchAddress("RockComp.second", RockComp_second, &b_RockComp_second);
   fChain->SetBranchAddress("MantleComp", &MantleComp_, &b_Header_MantleComp_);
   fChain->SetBranchAddress("MantleComp.first", MantleComp_first, &b_MantleComp_first);
   fChain->SetBranchAddress("MantleComp.second", MantleComp_second, &b_MantleComp_second);
   fChain->SetBranchAddress("CoreComp", &CoreComp_, &b_Header_CoreComp_);
   fChain->SetBranchAddress("CoreComp.first", CoreComp_first, &b_CoreComp_first);
   fChain->SetBranchAddress("CoreComp.second", CoreComp_second, &b_CoreComp_second);
   fChain->SetBranchAddress("RangeSW.fUniqueID", &RangeSW_fUniqueID, &b_Header_RangeSW_fUniqueID);
   fChain->SetBranchAddress("RangeSW.fBits", &RangeSW_fBits, &b_Header_RangeSW_fBits);
   fChain->SetBranchAddress("RangeSW.fName", &RangeSW_fName, &b_Header_RangeSW_fName);
   fChain->SetBranchAddress("RangeSW.fTitle", &RangeSW_fTitle, &b_Header_RangeSW_fTitle);
   fChain->SetBranchAddress("RangeSW.fLineColor", &RangeSW_fLineColor, &b_Header_RangeSW_fLineColor);
   fChain->SetBranchAddress("RangeSW.fLineStyle", &RangeSW_fLineStyle, &b_Header_RangeSW_fLineStyle);
   fChain->SetBranchAddress("RangeSW.fLineWidth", &RangeSW_fLineWidth, &b_Header_RangeSW_fLineWidth);
   fChain->SetBranchAddress("RangeSW.fFillColor", &RangeSW_fFillColor, &b_Header_RangeSW_fFillColor);
   fChain->SetBranchAddress("RangeSW.fFillStyle", &RangeSW_fFillStyle, &b_Header_RangeSW_fFillStyle);
   fChain->SetBranchAddress("RangeSW.fMarkerColor", &RangeSW_fMarkerColor, &b_Header_RangeSW_fMarkerColor);
   fChain->SetBranchAddress("RangeSW.fMarkerStyle", &RangeSW_fMarkerStyle, &b_Header_RangeSW_fMarkerStyle);
   fChain->SetBranchAddress("RangeSW.fMarkerSize", &RangeSW_fMarkerSize, &b_Header_RangeSW_fMarkerSize);
   fChain->SetBranchAddress("RangeSW.fNpoints", &RangeSW_fNpoints, &b_Header_RangeSW_fNpoints);
   fChain->SetBranchAddress("RangeSW.fX", RangeSW_fX, &b_RangeSW_fX);
   fChain->SetBranchAddress("RangeSW.fY", RangeSW_fY, &b_RangeSW_fY);
   fChain->SetBranchAddress("RangeSW.fMinimum", &RangeSW_fMinimum, &b_Header_RangeSW_fMinimum);
   fChain->SetBranchAddress("RangeSW.fMaximum", &RangeSW_fMaximum, &b_Header_RangeSW_fMaximum);
   fChain->SetBranchAddress("RangeSR.fUniqueID", &RangeSR_fUniqueID, &b_Header_RangeSR_fUniqueID);
   fChain->SetBranchAddress("RangeSR.fBits", &RangeSR_fBits, &b_Header_RangeSR_fBits);
   fChain->SetBranchAddress("RangeSR.fName", &RangeSR_fName, &b_Header_RangeSR_fName);
   fChain->SetBranchAddress("RangeSR.fTitle", &RangeSR_fTitle, &b_Header_RangeSR_fTitle);
   fChain->SetBranchAddress("RangeSR.fLineColor", &RangeSR_fLineColor, &b_Header_RangeSR_fLineColor);
   fChain->SetBranchAddress("RangeSR.fLineStyle", &RangeSR_fLineStyle, &b_Header_RangeSR_fLineStyle);
   fChain->SetBranchAddress("RangeSR.fLineWidth", &RangeSR_fLineWidth, &b_Header_RangeSR_fLineWidth);
   fChain->SetBranchAddress("RangeSR.fFillColor", &RangeSR_fFillColor, &b_Header_RangeSR_fFillColor);
   fChain->SetBranchAddress("RangeSR.fFillStyle", &RangeSR_fFillStyle, &b_Header_RangeSR_fFillStyle);
   fChain->SetBranchAddress("RangeSR.fMarkerColor", &RangeSR_fMarkerColor, &b_Header_RangeSR_fMarkerColor);
   fChain->SetBranchAddress("RangeSR.fMarkerStyle", &RangeSR_fMarkerStyle, &b_Header_RangeSR_fMarkerStyle);
   fChain->SetBranchAddress("RangeSR.fMarkerSize", &RangeSR_fMarkerSize, &b_Header_RangeSR_fMarkerSize);
   fChain->SetBranchAddress("RangeSR.fNpoints", &RangeSR_fNpoints, &b_Header_RangeSR_fNpoints);
   fChain->SetBranchAddress("RangeSR.fX", RangeSR_fX, &b_RangeSR_fX);
   fChain->SetBranchAddress("RangeSR.fY", RangeSR_fY, &b_RangeSR_fY);
   fChain->SetBranchAddress("RangeSR.fMinimum", &RangeSR_fMinimum, &b_Header_RangeSR_fMinimum);
   fChain->SetBranchAddress("RangeSR.fMaximum", &RangeSR_fMaximum, &b_Header_RangeSR_fMaximum);
   fChain->SetBranchAddress("detector", &detector, &b_Header_detector);
   fChain->SetBranchAddress("NuList", &NuList, &b_Header_NuList);
   fChain->SetBranchAddress("FluxFiles", &FluxFiles_, &b_Header_FluxFiles_);
   fChain->SetBranchAddress("FluxFiles.fUniqueID", FluxFiles_fUniqueID, &b_FluxFiles_fUniqueID);
   fChain->SetBranchAddress("FluxFiles.fBits", FluxFiles_fBits, &b_FluxFiles_fBits);
   fChain->SetBranchAddress("FluxFiles.NuType", FluxFiles_NuType, &b_FluxFiles_NuType);
   fChain->SetBranchAddress("FluxFiles.FluxSimul", FluxFiles_FluxSimul, &b_FluxFiles_FluxSimul);
   fChain->SetBranchAddress("FluxFiles.FileNames", FluxFiles_FileNames, &b_FluxFiles_FileNames);
   fChain->SetBranchAddress("PropMode", &PropMode, &b_Header_PropMode);
   fChain->SetBranchAddress("Drawing", &Drawing, &b_Header_Drawing);
   fChain->SetBranchAddress("gSeaGenVer", &gSeaGenVer, &b_Header_gSeaGenVer);
   fChain->SetBranchAddress("GenieVer", &GenieVer, &b_Header_GenieVer);
   fChain->SetBranchAddress("RunTime", &RunTime, &b_Header_RunTime);
   Notify();
}

Bool_t GSGHeaderReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GSGHeaderReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GSGHeaderReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GSGHeaderReader_cxx
