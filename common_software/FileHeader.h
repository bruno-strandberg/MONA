#ifndef FileHeader_h
#define FileHeader_h

#include "TNamed.h"
#include "TFile.h"
#include <map>
#include <vector>

class FileHeader;

/** A class that provides a minimalistic header functionality for ROOT files with easy access.
 *
 *  The class can be used to create 'headers' for any .root file. The idea is simple: it creates a
 *  directory named fHeaderDir (currently hardcoded to 'Header') into the output root file and writes
 *  TNamed objects into this directory, where the name of the TNamed object is in format
 *  ApplicationName__ParameterName__ParameterValue ('__' acts as delim and is configured through
 *  variable fDelim). In this way the header fields can be read without any extra software by just
 *  exploring the root file with TBrowser or in terminal.
 *
 *  Example usage:
 *  FileHeader a("myapp");                    # create header for application 'myapp'
 *  FileHeader.AddParameter("Rcan", "200");   # add parameter for can radius
 *  TFile fout("test.root","RECREATE");       # root file output
 *  //..some other stuff to be written to file..
 *  a.WriteHeader();                          # adds the header
 *  fout.Close();
 * 
 */
class FileHeader {

 public:
  FileHeader(TString application_name);
  ~FileHeader();

  void    AddParameter(TString parameter_name, TString parameter_value);
  TString GetParameter(TString application_name, TString parameter_name);

  void ReadHeader(TString filename);
  void WriteHeader(TFile *f);
  void AddToFile(TString filename);
  void Print();

 private:
  /// Name of the application where Header is initiated
  TString fAppName;
  /// Name of the dir in output root file where header is written (hardcoded to 'Header' in constructor)
  TString fHeaderDir;
  /// Delimiter (hardcoded to '__' in constructor)
  TString fDelim;
 //!< map of parameters for function GetParameter
  std::map< TString, std::map<TString, TString> > fPars;

  void CheckName(TString name);
  std::pair< TString, std::pair<TString, TString> > SplitDataStr(TString datastr);
};

#endif
