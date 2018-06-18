#ifndef SummaryParser_h
#define SummaryParser_h

#include <TChain.h>
#include <TFile.h>
#include "SummaryEvent.h"

#include <iostream>
using namespace std;

/**
   A small class to read data in SummaryEvent format.
 */
class SummaryParser {
  
 private:
  TChain         *fChain;   //!< pointer to the TChain where input files are attached
  SummaryEvent   *fEvt;     //!< pointer to the event structure
  
 public:
   SummaryParser(TString fname="");
   ~SummaryParser();

   TChain*        GetTree(); //!< tchain fetcher
   SummaryEvent*  GetEvt();  //!< event structure fetcher

};

#endif
