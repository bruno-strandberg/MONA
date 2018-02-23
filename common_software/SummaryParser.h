//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 22 09:13:22 2018 by ROOT version 6.06/04
// from TTree summary/ORCA MC summary tree
// found on file: ../data/ORCA_MC_summary_all.root
//////////////////////////////////////////////////////////

#ifndef SummaryParser_h
#define SummaryParser_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>
using namespace std;

// Header file for the classes stored in the TTree if any.

class SummaryParser {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        MC_runID;
   Double_t        MC_evtID;
   Double_t        MC_w2;
   Double_t        MC_erange_start;
   Double_t        MC_is_CC;
   Double_t        MC_is_neutrino;
   Double_t        MC_type;
   Double_t        MC_dir_x;
   Double_t        MC_dir_y;
   Double_t        MC_dir_z;
   Double_t        MC_pos_x;
   Double_t        MC_pos_y;
   Double_t        MC_pos_z;
   Double_t        MC_energy;
   Double_t        MC_bjorkeny;
   Double_t        gandalf_dir_x;
   Double_t        gandalf_dir_y;
   Double_t        gandalf_dir_z;
   Double_t        gandalf_pos_x;
   Double_t        gandalf_pos_y;
   Double_t        gandalf_pos_z;
   Double_t        gandalf_energy_nu;
   Double_t        shower_dir_x;
   Double_t        shower_dir_y;
   Double_t        shower_dir_z;
   Double_t        shower_pos_x;
   Double_t        shower_pos_y;
   Double_t        shower_pos_z;
   Double_t        shower_energy_nu;
   Double_t        shower_bjorkeny;
   Double_t        recolns_dir_x;
   Double_t        recolns_dir_y;
   Double_t        recolns_dir_z;
   Double_t        recolns_pos_x;
   Double_t        recolns_pos_y;
   Double_t        recolns_pos_z;
   Double_t        recolns_energy_nu;
   Double_t        recolns_bjorkeny;
   Double_t        PID_muon_probability;
   Double_t        PID_track_probability;

   // List of branches
   TBranch        *b_MC_runID;   //!
   TBranch        *b_MC_evtID;   //!
   TBranch        *b_MC_w2;   //!
   TBranch        *b_MC_erange_start;   //!
   TBranch        *b_MC_is_CC;   //!
   TBranch        *b_MC_is_neutrino;   //!
   TBranch        *b_MC_type;   //!
   TBranch        *b_MC_dir_x;   //!
   TBranch        *b_MC_dir_y;   //!
   TBranch        *b_MC_dir_z;   //!
   TBranch        *b_MC_pos_x;   //!
   TBranch        *b_MC_pos_y;   //!
   TBranch        *b_MC_pos_z;   //!
   TBranch        *b_MC_energy;   //!
   TBranch        *b_MC_bjorkeny;   //!
   TBranch        *b_gandalf_dir_x;   //!
   TBranch        *b_gandalf_dir_y;   //!
   TBranch        *b_gandalf_dir_z;   //!
   TBranch        *b_gandalf_pos_x;   //!
   TBranch        *b_gandalf_pos_y;   //!
   TBranch        *b_gandalf_pos_z;   //!
   TBranch        *b_gandalf_energy_nu;   //!
   TBranch        *b_shower_dir_x;   //!
   TBranch        *b_shower_dir_y;   //!
   TBranch        *b_shower_dir_z;   //!
   TBranch        *b_shower_pos_x;   //!
   TBranch        *b_shower_pos_y;   //!
   TBranch        *b_shower_pos_z;   //!
   TBranch        *b_shower_energy_nu;   //!
   TBranch        *b_shower_bjorkeny;   //!
   TBranch        *b_recolns_dir_x;   //!
   TBranch        *b_recolns_dir_y;   //!
   TBranch        *b_recolns_dir_z;   //!
   TBranch        *b_recolns_pos_x;   //!
   TBranch        *b_recolns_pos_y;   //!
   TBranch        *b_recolns_pos_z;   //!
   TBranch        *b_recolns_energy_nu;   //!
   TBranch        *b_recolns_bjorkeny;   //!
   TBranch        *b_PID_muon_probability;   //!
   TBranch        *b_PID_track_probability;   //!

   SummaryParser(TString fname);
   virtual ~SummaryParser();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   //------------------------------------------------------------
   //bstr: user defined methods, variables, etc
   //------------------------------------------------------------
};

#endif

#ifdef SummaryParser_cxx
SummaryParser::SummaryParser(TString fname) : fChain(0) 
{

  TFile *f = new TFile(fname, "READ");
  TTree *t = NULL;
  
  if ( f->IsOpen() ) { t = (TTree*)f->Get("summary"); }
  else { cout << "ERROR! SummaryParser::SummaryParser() coult not open file " << fname << endl; }

  if (t != NULL) Init(t);
  else { cout << "ERROR! SummaryParser::SummaryParser() could not find TTree 'summary' in "
	      << fname << endl; }
}

SummaryParser::~SummaryParser()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SummaryParser::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SummaryParser::LoadTree(Long64_t entry)
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

void SummaryParser::Init(TTree *tree)
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

   fChain->SetBranchAddress("MC_runID", &MC_runID, &b_MC_runID);
   fChain->SetBranchAddress("MC_evtID", &MC_evtID, &b_MC_evtID);
   fChain->SetBranchAddress("MC_w2", &MC_w2, &b_MC_w2);
   fChain->SetBranchAddress("MC_erange_start", &MC_erange_start, &b_MC_erange_start);
   fChain->SetBranchAddress("MC_is_CC", &MC_is_CC, &b_MC_is_CC);
   fChain->SetBranchAddress("MC_is_neutrino", &MC_is_neutrino, &b_MC_is_neutrino);
   fChain->SetBranchAddress("MC_type", &MC_type, &b_MC_type);
   fChain->SetBranchAddress("MC_dir_x", &MC_dir_x, &b_MC_dir_x);
   fChain->SetBranchAddress("MC_dir_y", &MC_dir_y, &b_MC_dir_y);
   fChain->SetBranchAddress("MC_dir_z", &MC_dir_z, &b_MC_dir_z);
   fChain->SetBranchAddress("MC_pos_x", &MC_pos_x, &b_MC_pos_x);
   fChain->SetBranchAddress("MC_pos_y", &MC_pos_y, &b_MC_pos_y);
   fChain->SetBranchAddress("MC_pos_z", &MC_pos_z, &b_MC_pos_z);
   fChain->SetBranchAddress("MC_energy", &MC_energy, &b_MC_energy);
   fChain->SetBranchAddress("MC_bjorkeny", &MC_bjorkeny, &b_MC_bjorkeny);
   fChain->SetBranchAddress("gandalf_dir_x", &gandalf_dir_x, &b_gandalf_dir_x);
   fChain->SetBranchAddress("gandalf_dir_y", &gandalf_dir_y, &b_gandalf_dir_y);
   fChain->SetBranchAddress("gandalf_dir_z", &gandalf_dir_z, &b_gandalf_dir_z);
   fChain->SetBranchAddress("gandalf_pos_x", &gandalf_pos_x, &b_gandalf_pos_x);
   fChain->SetBranchAddress("gandalf_pos_y", &gandalf_pos_y, &b_gandalf_pos_y);
   fChain->SetBranchAddress("gandalf_pos_z", &gandalf_pos_z, &b_gandalf_pos_z);
   fChain->SetBranchAddress("gandalf_energy_nu", &gandalf_energy_nu, &b_gandalf_energy_nu);
   fChain->SetBranchAddress("shower_dir_x", &shower_dir_x, &b_shower_dir_x);
   fChain->SetBranchAddress("shower_dir_y", &shower_dir_y, &b_shower_dir_y);
   fChain->SetBranchAddress("shower_dir_z", &shower_dir_z, &b_shower_dir_z);
   fChain->SetBranchAddress("shower_pos_x", &shower_pos_x, &b_shower_pos_x);
   fChain->SetBranchAddress("shower_pos_y", &shower_pos_y, &b_shower_pos_y);
   fChain->SetBranchAddress("shower_pos_z", &shower_pos_z, &b_shower_pos_z);
   fChain->SetBranchAddress("shower_energy_nu", &shower_energy_nu, &b_shower_energy_nu);
   fChain->SetBranchAddress("shower_bjorkeny", &shower_bjorkeny, &b_shower_bjorkeny);
   fChain->SetBranchAddress("recolns_dir_x", &recolns_dir_x, &b_recolns_dir_x);
   fChain->SetBranchAddress("recolns_dir_y", &recolns_dir_y, &b_recolns_dir_y);
   fChain->SetBranchAddress("recolns_dir_z", &recolns_dir_z, &b_recolns_dir_z);
   fChain->SetBranchAddress("recolns_pos_x", &recolns_pos_x, &b_recolns_pos_x);
   fChain->SetBranchAddress("recolns_pos_y", &recolns_pos_y, &b_recolns_pos_y);
   fChain->SetBranchAddress("recolns_pos_z", &recolns_pos_z, &b_recolns_pos_z);
   fChain->SetBranchAddress("recolns_energy_nu", &recolns_energy_nu, &b_recolns_energy_nu);
   fChain->SetBranchAddress("recolns_bjorkeny", &recolns_bjorkeny, &b_recolns_bjorkeny);
   fChain->SetBranchAddress("PID_muon_probability", &PID_muon_probability, &b_PID_muon_probability);
   fChain->SetBranchAddress("PID_track_probability", &PID_track_probability, &b_PID_track_probability);
   Notify();
}

Bool_t SummaryParser::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SummaryParser::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SummaryParser::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SummaryParser_cxx
