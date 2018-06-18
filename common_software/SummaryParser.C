#include "SummaryParser.h"
#include <stdexcept>

/**
   Constructor.
   \param fname     Name of the file to be parsed; wildcard accepted to add multiple files if read mode. In write mode the name of the file to be created where the data is written.
   \param ReadMode  If true the file opened for reading; if false the file is opened for writing
 */
SummaryParser::SummaryParser(TString fname, Bool_t ReadMode) {

  if (fname == "") {
    throw std::invalid_argument( "ERROR! SummaryParser::SummaryParser() empty filename as argument." );
  }
  
  fOut   = NULL;
  fEvt   = new SummaryEvent();
  fReadMode = ReadMode;

  if (fReadMode) {
    fChain = new TChain("summary");
    fChain->Add(fname);
    fChain->SetBranchAddress("SummaryEvent", &fEvt);  
  }
  else {
    fOut = new TFile(fname, "RECREATE");
    fChain = (TChain*)new TTree("summary", "ORCA MC summary tree");
    fChain->Branch("SummaryEvent", &fEvt, 2); //split level two means the tree is flattened
  }

}

/**
   Destructor.
 */
SummaryParser::~SummaryParser() {

  delete fEvt;
  if (fReadMode) {
    delete fChain;
  }
  else {
    if ( fOut->IsOpen() ) WriteAndClose(); //write tree and close file if not already done
    if (fOut) delete fOut; //fChain deleted when fOut is closed (I think...)
  }
  
}
