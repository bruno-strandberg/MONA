#include "PMNS_Base.h"
#include "PMNS_Deco.h"
#include "FitPDF.h"
#include "FitUtilWsyst.h"
#include "EffMass.h"
#include "DetResponse.h"

#include <iostream>
using namespace std;

/** class that does nothing with fProb */
class FitUtilNSI_1 : public FitUtilWsyst {

public:

  FitUtilNSI_1(Double_t op_time, TH3 *h_template, Double_t emin, Double_t emax, 
	     Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax, TString meff_file) : 
    FitUtilWsyst(op_time, h_template, emin, emax, ctmin, ctmax, bymin, bymax, meff_file) {

  }

  virtual ~FitUtilNSI_1() {};

};

//======================================================================================

/** class that creates a new fProb but does not overload ConfigOscProb */
class FitUtilNSI_2 : public FitUtilWsyst {

public:

  FitUtilNSI_2(Double_t op_time, TH3 *h_template, Double_t emin, Double_t emax, 
	     Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax, TString meff_file) : 
    FitUtilWsyst(op_time, h_template, emin, emax, ctmin, ctmax, bymin, bymax, meff_file) {

    if (fProb) delete fProb;
    fProb = new OscProb::PMNS_Deco;

  }

  virtual ~FitUtilNSI_2() {};

};

//======================================================================================

/** class that creates a new fProb and overloads ConfigOscProb */
class FitUtilNSI_3 : public FitUtilWsyst {

public:

  FitUtilNSI_3(Double_t op_time, TH3 *h_template, Double_t emin, Double_t emax, 
	     Double_t ctmin, Double_t ctmax, Double_t bymin, Double_t bymax, TString meff_file) : 
    FitUtilWsyst(op_time, h_template, emin, emax, ctmin, ctmax, bymin, bymax, meff_file) {

    if (fProb) delete fProb;
    fProb = new OscProb::PMNS_Deco;

  }

  virtual Bool_t ConfigOscProb(const proxymap_t& proxymap) { return kFALSE; }
  virtual ~FitUtilNSI_3() {};

};

//======================================================================================

/** This test checks that FitUtil behaves as expected when FitUtil::fProb is overloaded with a different OscProb calculator */
int main() {
  
  DetResponse dummy(DetResponse::track, "dummy");

  // init the utilities
  FitUtilNSI_1 NSI1(3, dummy.GetHist3D(), 3,60,-1,0, 0, 1, EffMass::DUMMYFILE);
  FitUtilNSI_2 NSI2(3, dummy.GetHist3D(), 3,60,-1,0, 0, 1, EffMass::DUMMYFILE);
  FitUtilNSI_3 NSI3(3, dummy.GetHist3D(), 3,60,-1,0, 0, 1, EffMass::DUMMYFILE);

  // init associated pdf's for access to proxymap
  FitPDF pdf1("pdf1","pdf1", &NSI1, &dummy);
  FitPDF pdf2("pdf2","pdf2", &NSI2, &dummy);
  FitPDF pdf3("pdf3","pdf3", &NSI3, &dummy);

  // create true bin object
  TrueB tb(1, 1, 0, 10, 5, 1);

  // check that NSI1 and NSI3 can be executed
  NSI1.TrueEvts( tb, pdf1.GetProxyMap() );
  NSI3.TrueEvts( tb, pdf3.GetProxyMap() );

  // check that NSI2 throw's a logic error
  Bool_t caught = kFALSE;

  try {
    NSI2.TrueEvts( tb, pdf2.GetProxyMap() );
  }
  catch (const std::logic_error& le) {
    caught = kTRUE;
  }
  
  if (!caught) {
    cout << "NOTICE TestProbOverload test failed" << endl;
    return 1;
  }
  else {
    cout << "NOTICE TestProbOverload test passed" << endl;
    return 0;
  }

}
