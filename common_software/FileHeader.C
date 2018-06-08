#include "FileHeader.h"
#include <iostream>
#include <stdexcept>

using namespace std;

//************************************************************************************

/**
 * Constructor
 * \param application_name Name of the application where the header is created.
 */
FileHeader::FileHeader(TString application_name) {

  fAppName   = application_name;
  fHeaderDir = "Header";
  fDelim     = "__";

  CheckName(application_name);

  fPars.insert( std::make_pair( application_name, std::map<TString, TString>() ) );

}

//************************************************************************************

/** Destructor */
FileHeader::~FileHeader() {}

//************************************************************************************

/** Function that prints out the header */
void FileHeader::Print() {
  
  for (auto &x: fPars) {
    cout << "Application: '" << x.first << "' parameters: " << endl;
    for (auto &y: x.second) {
      cout << "\t" << y.first << "\t" << y.second << endl;
    }
  }

}

//************************************************************************************

/** 
 * Function to add a parameter associated with the current application.
 *
 * \param parameter_name   Name of the parameter
 * \param parameter_value  Value of the parameter as a TString ( e.g. (TString)to_string(0.4) ).
*/
void FileHeader::AddParameter(TString parameter_name, TString parameter_value) {

  CheckName(parameter_name);
  CheckName(parameter_value);

  fPars[fAppName].insert( std::make_pair(parameter_name, parameter_value) );
  
}

//************************************************************************************

/**
 * Function to fetch the parameter from the header of a certain application.
 *
 * \param   application_name Name of the application the parameter is associated with
 * \param   parameter_name   Name of the parameter
 * \return  returns a TString with the parameter value, if found; otherwise empty string
 */
TString FileHeader::GetParameter(TString application_name, TString parameter_name) {

  TString par_value = "";

  if ( fPars.find(application_name) != fPars.end() ) {
    if ( fPars[application_name].find(parameter_name) != fPars[application_name].end() ) {
      par_value = fPars[application_name][parameter_name];
    }
  }

  if (par_value == "") {
    cout << "WARNING! FileHeader::GetParameter() could not find parameter " 
	 << parameter_name << " for application " << application_name << endl;
  }

  return par_value;
}

//************************************************************************************

/**
 * Function to write the header to an open root file (TFile).
 * Example: TFile fout("test.root","RECREATE"); fh.WriteHeader(&fout); fh.Close();
 * \param f Pointer to the opened root file where header should be written.
 */
void FileHeader::WriteHeader(TFile *f) {

  if (f == NULL) {
    throw std::invalid_argument( "ERROR! FileHeader::Write() null pointer to the file where header to be written." );
  }

  vector<TNamed> list;
  for (auto &app: fPars) {
    for (auto &par: app.second) {
      TString datastr = app.first + fDelim + par.first + fDelim + par.second;
      list.push_back( TNamed(datastr, (TString)"header field") );
    }
  }

  TDirectory *d = f->mkdir(fHeaderDir);
  d->cd();
  for (auto &l: list) l.Write();
  d->cd();

}

//************************************************************************************

/**
 * Function to add the header to an existing root file; if the file already contains a header, 
 * it is replaced.
 *
 * \param filename Name of the file where the header is written.
 */
void FileHeader::AddToFile(TString filename) {

  TFile f(filename, "update");

  // if previous header exists, delete it
  if ( f.Get(fHeaderDir) != NULL ) {
    f.Delete( fHeaderDir + TString(";*") );
  }

  // write this header to the file
  WriteHeader(&f);
  f.Close();
    
}

//************************************************************************************

/**
 * Function to read the header from another file to this header instance.
 * Example: FileHeader h("myapp"); h.ReadHeader("some_file.root")
 *
 * \param filename Name of the root file from where the header should be read.
 */
void FileHeader::ReadHeader(TString filename) {

  TFile f(filename, "READ");

  if ( !f.IsOpen() ) {
    throw std::invalid_argument("ERROR! FileHeader::ReadHeader() cannot open file " + (string)filename);
  }

  TDirectory *d = f.GetDirectory(fHeaderDir);
  if ( d == NULL ) {
    cout << "WARNING! FileHeader::ReadHeader() cannot find directory " << fHeaderDir << " in file " << filename << ", header reading failed." << endl;
    return;
  }

  TList *l = d->GetListOfKeys();

  for (auto x: *l) {
    TNamed *hf = (TNamed*)d->Get( x->GetName() );
    auto data = SplitDataStr( hf->GetName() );
    fPars[data.first].insert( std::make_pair(data.second.first, data.second.second) );    
  }
  
}

//************************************************************************************

/**
 *  Private function to check that the application name, parameter name and parameter value 
 *  string conform, throws an exception if they dont.
 *
 *  \param name Name of the application, parameter name or parameter value.
 */
void FileHeader::CheckName(TString name) {

  if (name == "") {
    throw std::invalid_argument( "ERROR! FileHeader::CheckName() empty string not allowed" );
  }

  if ( name.Contains(fDelim) ) {
    throw std::invalid_argument( "ERROR! FileHeader::CheckName() character combination '" + (string)fDelim + "' is reserved as a delimiter, change " + (string)name );
  }

}

//************************************************************************************

/**
 * Private function to split a string read from a header file to fields.
 
 * \param datastr Assuming deliminator '__', string in format Application__ParName__ParValue
 * \return a pair, where first element is Application, second element is a pair (ParName, ParValue)
 */
std::pair< TString, std::pair<TString, TString> > FileHeader::SplitDataStr(TString datastr) {

  TString remainder = datastr;
  vector<TString> parts;
  
  while ( remainder.Contains(fDelim) ) {    
    Int_t idx = remainder.Index(fDelim);
    parts.push_back( (TString)remainder(0, idx) );
    remainder = (TString)remainder( idx + fDelim.Length(), remainder.Length() );
  }

  parts.push_back(remainder);

  if (parts.size() != 3) {
    for (auto &s: parts) cout << "\t" << s << endl;
    throw std::invalid_argument("ERROR! FileHeader::SplitDataStr() expected 3 parts in " + (string)datastr + ", but got the " + to_string(parts.size()) + " printed above");
  }

  return std::make_pair( parts[0], std::make_pair(parts[1], parts[2]) );

}
