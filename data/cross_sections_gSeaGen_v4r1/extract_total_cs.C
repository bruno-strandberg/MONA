/*
  Extract the TGraphs containing the total interaction cross-section.
*/

void copyTGraph( TFile* fin, TFile* fout, string dir, string gname ) {
  string fullname = dir + "/" + gname ;
  TGraph* g = (TGraph*) fin->Get( fullname.c_str() ) ;

  if ( g== NULL ) {
    cerr << "FATAL ERROR: could not get TGraph '" 
	 << gname
	 << "' from directory '"
	 << dir
	 << "'." << endl ;
    exit(1) ;
  }

  fout->cd( dir.c_str() ) ;
  g->Write() ;
}

void extract_total_cs() {

  const string ifname = "xsec_full.root" ;
  const string ofname = "xsec.root" ;

  TFile* fin  = new TFile( ifname.c_str(), "read" ) ;
  TFile* fout = new TFile( ofname.c_str(), "recreate" ) ;

  const int nflavours = 3 ;
  string flavours[nflavours] = { "e", "mu", "tau" } ;

  const int ntargets = 3 ;
  string targets[ntargets] = { "n", "H1", "O16" } ;


  // loop over flavours
  for( int iflavour=0; iflavour!=nflavours; ++iflavour ) {
    string flavour = flavours[iflavour] ;
    string dir1 = "nu_" + flavour ;

    // loop over neutrino/antineutrino
    for( int ip=0 ; ip<2; ++ip ) {
      if( ip == 1 ) dir1 += "_bar" ;

      // loop over targets
      for( int itarget=0; itarget!=ntargets; ++itarget ) {
	string target = targets[itarget] ;
        string dir = dir1 + "_" + target ;

	fout->mkdir( dir.c_str() ) ;

	copyTGraph( fin, fout, dir, "tot_nc" ) ;
	copyTGraph( fin, fout, dir, "tot_cc" ) ;
      
      }

    }
  }

  fin->Close() ;
  fout->Close() ;

  cout << "Output written to '" << ofname << "'." << endl ;

}
