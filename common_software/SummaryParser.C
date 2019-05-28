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
  fHead  = new FileHeader("SummaryParser");
  fReadMode = ReadMode;

  TString SumEvtV = "SummaryEventVersion";                       // parameter name for the header
  Int_t   SEV_lib = ( (TClass*)fEvt->IsA() )->GetClassVersion(); // summary event version in the library
  
  if (fReadMode) {

    // get the header
    fHead->ReadHeader(fname);
    TString sevstr_file = fHead->GetParameter(SumEvtV);
    Int_t   SEV_file    = 0;

    // if it cannot find the header throw an error
    if ( sevstr_file != "" ) {
      SEV_file = stoi( (string)sevstr_file );
    }
    else {
      throw std::invalid_argument("ERROR! SummaryParser::SummaryParser() the input file " + (string)fname + " does not contain a header with the SummaryEvent version info.");
    }

    // if SummaryEvent versions differ throw an error
    if ( SEV_file != SEV_lib ) {
      throw std::invalid_argument("ERROR! SummaryParser::SummaryParser() the input file " + (string)fname + " contains SummaryEvent's of version " + to_string(SEV_file) + ", but the library has version " + to_string(SEV_lib) + ". Either check out an older version of MONA or create a new summary file with the tools in /apps/data_sorting." );
    }
    
    fChain = new TChain("summary");
    fChain->Add(fname);
    fChain->SetBranchAddress("SummaryEvent", &fEvt);  
  }
  else {
    fOut = new TFile(fname, "RECREATE");
    fChain = (TChain*)new TTree("summary", "ORCA MC summary tree");
    fChain->Branch("SummaryEvent", &fEvt, 2); //split level two means the tree is flattened
    fHead->AddParameter(SumEvtV, (TString)to_string(SEV_lib) );
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

  delete fHead;
  
}
