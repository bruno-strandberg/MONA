#ifndef SummaryParser_h
#define SummaryParser_h

#include <TChain.h>
#include <TFile.h>
#include "SummaryEvent.h"

#include <iostream>
using namespace std;

/**
   A small class to read/write data in SummaryEvent format.
 */
class SummaryParser {
  
 private:
  TChain         *fChain;     //!< pointer to the TChain where input files are attached
  SummaryEvent   *fEvt;       //!< pointer to the event structure
  TFile          *fOut;       //!< pointer to output file in write mode
  Bool_t          fReadMode;  //!< flag to indicate read/write mode
  
 public:
  SummaryParser(TString fname, Bool_t ReadMode=kTRUE);
   ~SummaryParser();

   /// write mode - write the filled tree to fOut and close file
   void           WriteAndClose() { fOut->cd(); fChain->Write(); fOut->Close(); };

   // getters

   /** Function to fetch the pointer to the `fChain` class member.
       \return Pointer to the `fChain` class member.
    */
   TChain*        GetTree()        { return fChain; }

   /** Function to fetch the pointer to the `fEvt` class member.

       NB! When using this, first do `GetTree()->GetEntry(i)` to read the event `i` to `fEvt`.

       \return Pointer to the `fEvt` class member.
      
    */
   SummaryEvent*  GetEvt()         { return fEvt;   }      //!< event structure fetcher

   /** Function to fetch the pointer to the `fEvt` class member, which also calls `GetEntry(i)` on `fChain`.
       \param i     Event number
       \return      Pointer to the `fEvt` class member corresponding to the event number `i`.
    */
   SummaryEvent*  GetEvt(Int_t i)  { if (fChain) fChain->GetEntry(i); return fEvt; }

};

#endif
