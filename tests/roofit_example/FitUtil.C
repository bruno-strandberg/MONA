#include "FitUtil.h"
#include <iostream>

using namespace std;

/** Constructor*/
FitUtil::FitUtil() {

  Ereco  = new RooRealVar("Ereco"  ,"Ereco" , 1, 0, 20);
  Ctreco = new RooRealVar("Ctreco" ,"Ctreco", 1, 0, 20);
  a      = new RooRealVar("a","a", 1, -10, 10);
  b      = new RooRealVar("b","b", 1, -10, 10);

  // all variables need to enter this list!
  fParSet.add( RooArgSet(*Ereco, *Ctreco, *a, *b) );

  // only observables need to enter this list!
  fObsList.add( RooArgList( *Ereco, *Ctreco ) );

}

//***************************************************************************

/** Destructor.*/
FitUtil::~FitUtil() {

  RooArgList l( fParSet );

  for (Int_t i = 0; i < l.getSize(); i++) {
    RooRealVar *v = (RooRealVar*)l.at(i);
    if (v) delete v;
  }

}

//***************************************************************************

/** Function that returns a value depending on the variables
    \parameter Ereco  Energy observable
    \parameter ctreco cos-theta observable
    \parameter a      model parameter a
    \parameter b      model parameter b
    \return           Model value
*/
Double_t FitUtil::GetValue(Double_t Ereco, Double_t ctreco, Double_t a, Double_t b) {
  return (a+b) * Ereco + (a-b) * ctreco;
}

//***************************************************************************

/** This function is called in `FitPDF::evaluate()`
    
    \parameter parmap   Map of parameter names and `RooRealProxy`'s
    \return             Model value

*/
Double_t FitUtil::GetValue(const std::map<TString, RooRealProxy*> &parmap) {

  Double_t energy = *(parmap.at( (TString)Ereco->GetName()  ) );
  Double_t ct     = *(parmap.at( (TString)Ctreco->GetName() ) );
  Double_t a_     = *(parmap.at( (TString)a->GetName() ) );
  Double_t b_     = *(parmap.at( (TString)b->GetName() ) );

  return GetValue( energy, ct, a_, b_ );

}
