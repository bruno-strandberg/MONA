#include "FileHeader.h"
#include "TSystem.h"
#include <iostream>
#include "TFile.h"

using namespace std;

void fileheader() {

  gSystem->Load("$NMHDIR/common_software/libnmhsoft.so");

  // create a header and add some parameters to it, print it out
  FileHeader a("fileheader");
  a.AddParameter( "sinsq_12", (TString)to_string(0.297) );
  a.AddParameter( "sinsq_23", (TString)to_string(45) );
  cout << "==================Printing header a=========================" << endl;
  a.Print();
  cout << "==================Finished==================================" << endl;

  // write the header to file
  TFile fh("fileheader.root","RECREATE");
  a.WriteHeader(&fh);
  fh.Close();

  // create another header, add some parameters, print
  FileHeader b("secondapp");
  b.AddParameter("Rcan", (TString)to_string(190) );
  b.AddParameter("Zmin", (TString)to_string(0) );
  b.AddParameter("Zmax", (TString)to_string(210) );
  cout << "==================Printing header b before read-in==========" << endl;
  b.Print();
  cout << "==================Finished==================================" << endl;
  cout << b.GetParameter("secondapp", "Rcan") << endl;

  // additionally read the header info from the previous file, print again
  b.ReadHeader("fileheader.root");
  cout << "==================Printing header b after read-in===========" << endl;
  b.Print();
  cout << "==================Finished==================================" << endl;
  cout << b.GetParameter("fileheader", "sinsq_12") << endl;

  // recreate the file and do one final read-in, read-out
  TFile fh2("fileheader.root","RECREATE");
  b.WriteHeader(&fh2);
  fh2.Close();

  FileHeader c("fileheader");
  c.ReadHeader("fileheader.root");
  cout << "==================Printing header c after===================" << endl;
  c.Print();
  cout << "==================Finished==================================" << endl;

}
