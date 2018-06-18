#include "SummaryParser.h"

/**
   Constructor.
   \param fname   Name of the file to be parsed.
 */
SummaryParser::SummaryParser(TString fname) {

  fChain = new TChain("summary", "ORCA MC summary tree");
  fEvt   = new SummaryEvent();
  fChain->SetBranchAddress("SummaryEvent", &fEvt);  
  
  if (fname != "") fChain->Add(fname);

}

/**
   Destructor.
 */
SummaryParser::~SummaryParser() {

  delete fEvt;
  delete fChain;

}
