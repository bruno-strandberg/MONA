using namespace std;

void fileheader() {

  // create a header and add some parameters to it, print it out
  FileHeader a("fileheader");
  a.AddParameter( "sinsq_12", (TString)to_string(0.297) );
  a.AddParameter( "sinsq_23", (TString)to_string(45) );
  a.AddParameter( "filename", "some/dir_here/structure.root" );
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

  // additionally read the header info from the previous file, print again
  b.ReadHeader("fileheader.root");
  cout << "==================Printing header b after read-in===========" << endl;
  b.Print();
  cout << "==================Finished==================================" << endl;

  cout << "==================Testing parameter finding=================" << endl;
  cout << "Sinsq_12 through get parameter            : " << b.GetParameter("sinsq_12","fileheader.root","fileheader") << endl;
  cout << "Sinsq_12 through get parameter (name only): " << b.GetParameter("sinsq_12") << endl;
  cout << "Sinsq_12 through find parameter           : " << b.FindParValues("sinsq_12").front().back() << endl;

  // recreate the file and do one more read-in, read-out
  TFile fh2("fileheader2.root","RECREATE");
  b.WriteHeader(&fh2);
  fh2.Close();

  FileHeader c("readin");
  c.ReadHeader("fileheader2.root");
  cout << "==================Printing header c after read-in===========" << endl;
  c.Print();
  cout << "==================Finished==================================" << endl;

  FileHeader d("updated");
  d.AddParameter("NewPar", "69");
  d.AddToFile("fileheader2.root", false);

  FileHeader e("readupdated");
  e.ReadHeader("fileheader2.root");
  cout << "==================Printing replaced header==================" << endl;
  e.Print();
  cout << "==================Finished==================================" << endl;  

  system("rm fileheader.root fileheader2.root");
}
