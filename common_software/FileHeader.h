#ifndef FileHeader_h
#define FileHeader_h

#include "TNamed.h"
#include "TFile.h"
#include <map>
#include <vector>

class FileHeader;

/** 
    A class that provides a minimalistic header functionality for ROOT files.
 
    The header is a `TTree` with character arrays as entries. The individual lines can be explored on the command line by e.g. `root myfile.root; HeaderTree->Show(i)`, where i = 0,1,...
    
    The header data is stored in a map with structure [outputname][application][parameter][value]. Outputname is the name of the `.root` output where the header originates(!) from, application is the name of the application or script where the header was created, parameter is the name of a parameter and value is the value of the parameter.
 
    With this structure, the header information of each output file can be preserved when:
    1) One uses a chain of applications
    2) One merges (e.g. with hadd) multiple outputs created with the same application.
    
    Example case 1: create a header into file fout.root which stores the can radius
    \code{.cpp}
    FileHeader a("myapp"); 
    a.AddParameter("Rcan","190"); 
    TFile f("fout.root","RECREATE"); 
    a.WriteHeader(&f); 
    f.Close()
    \endcode
    In this case the map will have a field fPars['fout.root']['myapp']['Rcan'] = '190'
 
    Example case 2: create a header, read header data from another file, write the header to another file
    \code{.cpp}
    FileHeader b("secondapp"); 
    b.AddParameter("mu_cut","0.05"); 
    b.ReadHeader("fout.root"); 
    TFile f("fnew.root","RECREATE"); 
    b.WriteHeader(&f),
    f.Close();
    \endcode
    In this case the map will have a field fPars['fout.root']['myapp']['Rcan'] = '190' and fPars['fnew.root']['secondapp']['mu_cut'] = '0.05'. Note that although header `b` is written to file `fnew.root`, the header will also contain data for output name `fout.root`. This allows to trace parameters through chains of applications.
   
    To print the header in the root command line, do
    \code{.cpp}
    root;
    FileHeader a;
    a.ReadHeader('myfile.root')
    a.Print()
    \endcode
*/
class FileHeader {

 public:
  FileHeader(TString application_name);
  ~FileHeader();

  void    AddParameter(TString parameter_name, TString parameter_value);
  void    AddParameter(FileHeader &h, TString par_name);
  TString GetParameter(TString parameter_name, TString output_name="", TString application_name="");
  std::vector< std::vector<TString> > FindParValues(TString parameter_name);

  void ReadHeader(TString filename);
  void WriteHeader(TFile *f);
  void AddToFile(TString filename, Bool_t overwrite=kFALSE);
  void Print();

 private:
  /// Name of the application where Header is initiated
  TString fAppName;
  /// Name of the tree where header is written (hardcoded to 'Header' in constructor)
  TString fHeaderTree;
  /// Delimiter (hardcoded to '__' in constructor)
  TString fDelim;
  /// Max line length (hardcoded to 1024 characters in constructor)
  UInt_t  fLineLength;
 //!< map of parameters for function GetParameter with structure fPars[outputname][application][parameter] = value
  std::map< TString, std::map< TString, std::map<TString, TString> > > fPars;

  void CheckName(TString name);
  std::vector<TString> SplitDataStr(TString datastr);
  void AddToMap(TString out_name, TString app_name, TString par_name="", TString par_value="");
};

#endif
