#include "FileHeader.h"
#include "TTree.h"
#include <iostream>
#include <stdexcept>

using namespace std;

//************************************************************************************

/**
 * Constructor
 * \param application_name Name of the application where the header is created.
 */
FileHeader::FileHeader(TString application_name) {

  fAppName    = application_name;
  fHeaderTree = "Header";
  fDelim      = "__";
  fLineLength = 1024;

  CheckName(application_name);

  // "this" is a placeholder for the filename where this header will be written to; it is replaced in WriteHeader()
  AddToMap("this", application_name);

}

//************************************************************************************

/** Destructor */
FileHeader::~FileHeader() {}

//************************************************************************************

/** Function that prints out the header */
void FileHeader::Print() {
  
  for (auto &x: fPars) {
    cout << "Output: " << x.first << endl;
    for (auto &y: x.second) {
      cout << "\t" << "Application: " << y.first << endl;
      cout << "\t\t" << "Parameters: " << endl;
      for (auto &z: y.second) {
	cout << "\t\t" << z.first << "\t" << z.second << endl;
      }
    }
  }

}

//************************************************************************************

/** 
 * Function to add a parameter associated with the current application and output
 *
 * \param parameter_name   Name of the parameter
 * \param parameter_value  Value of the parameter as a TString ( e.g. (TString)to_string(0.4) ).
*/
void FileHeader::AddParameter(TString parameter_name, TString parameter_value) {

  CheckName(parameter_name);
  CheckName(parameter_value);

  fPars["this"][fAppName].insert( std::make_pair(parameter_name, parameter_value) );
  
}

//************************************************************************************

/** 
 * Function to include parameter from another FileHeader in this header.
 *
 * This function is useful when one wishes to manually add a parameter from another header
 * instead of reading in the full header to this instance with ReadHeader(). If the parameter
 * name is present several times, the first found value is added.
 *
 * \param h          Reference to another FileHeader instance
 * \param par_name   Name of the parameter to be included
*/
void FileHeader::AddParameter(FileHeader &h, TString par_name) {

  CheckName(par_name);
  auto par = h.FindParValues(par_name);
  if ( par.size() > 0 ) {
    AddToMap("this", par.front()[1], par.front()[2], par.front()[3]);
  }

}

//************************************************************************************

/**
 * Function to fetch the parameter from the header of a certain output and application.
 *
 * If output_name and application_name are not specified, it will look through the map
 * and return the value of the first parameter with this name.
 *
 * \param   parameter_name   Name of the parameter
 * \param   output_name      Name of the output file the parameter is associated with
 * \param   application_name Name of the application the parameter is associated with
 * \return  returns a TString with the parameter value, if found; otherwise empty string
 */
TString FileHeader::GetParameter(TString parameter_name, TString output_name, TString application_name) {

  if (parameter_name == "") {
    throw std::invalid_argument("ERROR! FileHeader::GetParameter() empty string provided as parameter name");
  }
  
  TString par_value = "";

  if (output_name != "" && application_name != "") {

    if ( fPars.find(output_name) != fPars.end() ) {
      if ( fPars[output_name].find(application_name) != fPars[output_name].end() ) {
	if ( fPars[output_name][application_name].find(parameter_name) != fPars[output_name][application_name].end() ) {
	  par_value = fPars[output_name][application_name][parameter_name];
	}
      }
    }

  }
  else  {

    vector< vector<TString> > pars = FindParValues(parameter_name);
    if (pars.size() > 0) { par_value = pars.front().back(); }

  }
  
  if (par_value == "") {

    TString wrnmsg = "WARNING! FileHeader::GetParameter() could not find parameter " + parameter_name;
    if (output_name != "") wrnmsg += " for output " + output_name;
    if (application_name != "") wrnmsg += " for application " + application_name;

    cout << wrnmsg << endl;
    
  }

  return par_value;
}

//************************************************************************************

/** 
 * Function to look up parameters in a map by name only.
 * \param parameter_name Name of the parameter
 * \return A vector of vectors; each vector contains 4 TStrings: outname, appname, parname, parvalue
*/
vector< vector<TString> > FileHeader::FindParValues(TString parameter_name) {

  vector< vector<TString> > pars_vec;

  for (auto &outs: fPars) {
    for (auto &apps: outs.second) {
      for (auto &pars: apps.second) {
	if ( pars.first.Contains(parameter_name) ) {
	  vector<TString> par = { outs.first, apps.first, pars.first, pars.second };
	  pars_vec.push_back(par);
	}
      }
    }
  }

  return pars_vec;

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

  CheckName( f->GetName() );
  
  TTree *t = new TTree(fHeaderTree,fHeaderTree);
  Char_t line[fLineLength];
  t->Branch("HeaderLine", line, "HeaderLine/C");

  for (auto &out: fPars) {
    for (auto &app: out.second) {
      for (auto &par: app.second) {

	TString outname = f->GetName();
	TString appname = app.first;
	TString parname = par.first;
	TString valname = par.second;

	TString out_app_par_str = outname + fDelim + appname + fDelim + parname + fDelim + valname;

	sprintf(line, "%s", (const char*)out_app_par_str);
	t->Fill();
      }
    }
  }

  t->Write();

}

//************************************************************************************

/**
 * Function to add the header to an existing root file.
 *
 * \param filename  Name of the file where the header is written.
 * \param overwrite If set to true, the header is overwritten, otherwise data is added to the header.
 */
void FileHeader::AddToFile(TString filename, Bool_t overwrite) {

  if (!overwrite) ReadHeader(filename);

  TFile f(filename, "update");

  if ( !f.IsOpen() ) {
    throw std::invalid_argument("ERROR! FileHeader::AddToFile() cannot open file " + (string)filename);
  }

  // if previous header exists, delete it
  if ( f.Get(fHeaderTree) != NULL ) {
    gDirectory->Delete( fHeaderTree + TString(";*") );
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

  TTree *t = (TTree*)f.Get(fHeaderTree);
  if ( t == NULL ) {
    cout << "WARNING! FileHeader::ReadHeader() cannot find directory " << fHeaderTree << " in file " << filename << ", header reading failed." << endl;
    return;
  }

  Char_t line[fLineLength];
  t->SetBranchAddress("HeaderLine", line);

  for (Int_t i = 0; i < t->GetEntries(); i++) {
    t->GetEntry(i);
    vector<TString> data = SplitDataStr( (TString)line );
    AddToMap( data[0], data[1], data[2], data[3] );
  }

  f.Close();
  
}

//************************************************************************************

/**
 *  Private function to check that the output name,  application name, parameter name 
 *  and parameter value string conform, throws an exception if they dont.
 *
 *  \param name Name of the output, application, parameter name or parameter value.
 */
void FileHeader::CheckName(TString name) {

  if (name == "") {
    throw std::invalid_argument( "ERROR! FileHeader::CheckName() empty string not allowed" );
  }

  if ( name.Contains(fDelim) ) {
    throw std::invalid_argument( "ERROR! FileHeader::CheckName() character combination '" + (string)fDelim + "' is reserved as a delimiter, change " + (string)name );
  }

  if ( name.Length() > (Int_t)fLineLength) {
    throw std::invalid_argument( "ERROR! FileHeader::CheckName() name '" + (string)name + "' is longer than the maximum line length of " + to_string(fLineLength) );
  }

}

//************************************************************************************

/**
 * Private function to split a string read from a header file to fields.
 
 * \param datastr Assuming deliminator '__', string in format Output__Application__ParName__ParValue
 * \return a vector, where element [0] is output, [1] is Application, [2] is ParName, [3] is ParValue
 */
vector<TString> FileHeader::SplitDataStr(TString datastr) {

  TString remainder = datastr;
  vector<TString> parts;
  
  while ( remainder.Contains(fDelim) ) {    
    Int_t idx = remainder.Index(fDelim);
    parts.push_back( (TString)remainder(0, idx) );
    remainder = (TString)remainder( idx + fDelim.Length(), remainder.Length() );
  }

  parts.push_back(remainder);

  if (parts.size() != 4) {
    for (auto &s: parts) cout << "\t" << s << endl;
    throw std::invalid_argument("ERROR! FileHeader::SplitDataStr() expected 4 parts in " + (string)datastr + ", but got the " + to_string(parts.size()) + " printed above");
  }

  return parts;

}

//************************************************************************************

/** Private function to add elements to the map.
 *
 * \param out_name  Name of the output the parameter is associated with
 * \param app_name  Name of the application the parameter is associated with
 * \param par_name  Name of the parametre
 * \param par_value Value of the parameter as TString
 */
void FileHeader::AddToMap(TString out_name, TString app_name, TString par_name, TString par_value) {

  // check if fPars already has an entry for this out_name; if not, create it
  if ( fPars.find(out_name) == fPars.end() ) {
    
    // create a map application name <--> parameter map
    std::map<TString, std::map<TString, TString> > app_par_map;

    // create a pair output name <--> application parameter map
    fPars.insert( std::make_pair( out_name, app_par_map ) );

  }

  // insert will add an empty map for app_name; insert does not overwrite if map already exists
  fPars[out_name].insert( std::make_pair( app_name, std::map<TString, TString>() )  );

  // add the par_name, par_value pair to map; insert does not overwrite if map already exists
  if (par_name != "" && par_value != "") {
    fPars[out_name][app_name].insert( std::make_pair(par_name, par_value) );
  }
  else if ( (par_name == "" && par_value != "") || (par_name != "" && par_value == "") ) {
    throw std::invalid_argument( "ERROR! FileHeader::AddToMap() both par_name and par_value need to be set or need to be empty." );    
  }

}
