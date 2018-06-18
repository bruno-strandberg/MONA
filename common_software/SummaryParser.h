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
   TChain*        GetTree() { return fChain; }      //!< tchain fetcher
   SummaryEvent*  GetEvt()  { return fEvt;   }      //!< event structure fetcher

};

#endif
