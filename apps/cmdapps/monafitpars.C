
// MONA headers
#include "DetResponse.h"
#include "FitUtil.h"
#include "FitUtilWsyst.h"

//jpp headers
#include "JTools/JRange.hh"
#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

#include <stdexcept>

/** This application prints out the observable names and fit parameter names in the specified fit utility class in `MONA/fitter_software/` */
int main(const int argc, const char **argv) {

  TString FitUtil_str      = "FitUtil";
  TString FitUtilWsyst_str = "FitUtilWsyst";
  
  vector<TString> fuclasses = { FitUtil_str, FitUtilWsyst_str  };
  TString fitutility;

  try {

    JParser<> zap("This application prints out the observable names and fit parameter names in the specified fit utility class in MONA/fitter_software/");

    zap['f'] = make_field(fitutility, "Fit Utility class") = fuclasses;
    if ( zap.read(argc, argv) != 0 ) return 1;
    
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  DetResponse resp(DetResponse::track, "dummy");

  FitUtil *fu;

  if (fitutility == FitUtil_str) {
    fu = new FitUtil( 1, resp.GetHist3D(), 5, 50, -1, 0, 0, 1, EffMass::DUMMYFILE );
  }
  else if (fitutility == FitUtilWsyst_str) {
    fu = new FitUtilWsyst( 1, resp.GetHist3D(), 5, 50, -1, 0, 0, 1, EffMass::DUMMYFILE );
  }
  else {
    throw std::invalid_argument("ERROR! monafitpars: unknown fit utility class specified");
  }

  RooArgList pars( fu->GetSet() );
  TIterator *it = pars.createIterator();
  RooRealVar *var;

  vector< TString > parameters;
  vector< TString > observables;
  
  while ( ( var = (RooRealVar*)it->Next() ) ) {
    if ( fu->GetObs().find( var->GetName() ) ) {
      observables.push_back( (TString)var->GetName() );
    }
    else {
      parameters.push_back( (TString)var->GetName() );
    }
  }

  cout << "============================================================================" << endl;
  cout << "Observables (by name) in class " << fitutility << ":" << endl;
  for (auto o: observables) cout << o << endl;
  cout << "============================================================================" << endl;
  cout << "Fit parameters (by name) in class " << fitutility << ":" << endl;
  for (auto p: parameters) cout << p << endl;
  cout << "============================================================================" << endl;
  
}
