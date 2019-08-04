#include "EventFilter.h"
#include<iostream>
#include<stdexcept>
using namespace std;

/**
   Constructor.

   \param reco_type Type of reconstruction variables written to observables fEnergy, fDir, fPos, fBy, e.g. EventFilter::track.
 */
EventFilter::EventFilter(reco reco_type) {
  fRecoType = reco_type;

  // the pointers to custom observable functions are initialized to NULL and need to be configured
  // by the user, if "custom" reco type is used.
  GetCustomEnergy = NULL;
  GetCustomDir    = NULL;
  GetCustomPos    = NULL;
  GetCustomBy     = NULL;

  if (fRecoType == customreco) {
    std::cout << "NOTICE EventFilter::EventFilter() initialised with a customised reco type, call SetObsFuncPtrs() to configure the observables to be used." << std::endl;
  }

}

//***************************************************************************************

/** Destructor */
EventFilter::~EventFilter() {}

//***************************************************************************************

/**
   Function to add an 'and' cut or an 'or' cut to this event filter.

   See the class description for example usage.
   
   \param getter_ptr Pointer to a getter function of SummaryEvent class
   \param comp_ptr   Pointer to a comparison function.
   \param value      Value that the getter return is compared to
   \param AndCut     true - treat this cut as 'and' cut; false - treat this cut as 'or' cut
 */
void EventFilter::AddCut( std::function<Double_t(SummaryEvent&)> getter_ptr,
			  std::function<bool(double, double)> comp_ptr, 
			  Double_t value, Bool_t AndCut) {
  cutobj this_cut(getter_ptr, value, comp_ptr);
  if (AndCut) fAndCuts.push_back(this_cut);
  else fOrCuts.push_back(this_cut);
}

//***************************************************************************************

/**
   Function to add an 'and' cut or an 'or' cut to this event filter in python.

   Example usage:
   \code{.cpp}
   EventFilter f;
   f.AddCut("SummaryEvent.Get_track_ql0", ">", 0.5, kTRUE);
   \endcode
   
   The first string specifies a getter function name of `SummaryEvent` class, the second string specifies a comparison operator. When using MONA with PyROOT, only this function will work for adding cuts.

   \param getter_func_name  `SummaryEvent` function name, e.g. `SummaryEvent.Get_track_ql0`
   \param comp_op           Comparison operator string, supported types are ">", ">=", "==", "<=", "<"
   \param value             Value that the getter return is compared to
   \param AndCut            true - treat this cut as 'and' cut; false - treat this cut as 'or' cut
 */
void EventFilter::AddCut(TString getter_func_name, TString comp_op, Double_t value, Bool_t AndCut) {

  if ( SUMMARYUTILS::fFuncMap.find(getter_func_name) == SUMMARYUTILS::fFuncMap.end() ) {
    for (auto f: SUMMARYUTILS::fFuncMap) cout << f.first << endl;
    throw std::invalid_argument("ERROR! EventFilter::AddCut(TString ...) unrecognized function name " + (string)getter_func_name + " . Recognized function names listed above.");
  }

  if ( fOpMap.find(comp_op) == fOpMap.end() ) {
    for (auto o: fOpMap) cout << o.first << endl;
    throw std::invalid_argument("ERROR! EventFilter::AddCut(TString ...) unrecognized operator " + (string)comp_op + " . Recognized operators listed above.");
  }

  cutobj this_cut( SUMMARYUTILS::fFuncMap[getter_func_name], value, fOpMap[comp_op]);
  if (AndCut) fAndCuts.push_back(this_cut);
  else fOrCuts.push_back(this_cut);
}

//***************************************************************************************

/** 
    Function to determine whether the event passes the cuts defined for this event filter.

    The function performes an 'and' between all cuts added through AddAndCut and an 'or'
    between all cuts added through AddOrCut. It returns PassesAddCuts && PassesOrCuts.

    \param evt Pointer to a SummaryEvent
    \return    True if passes all cuts, false otherwise.
*/
Bool_t EventFilter::PassesCuts(SummaryEvent *evt) {

  /* As a comment: cut.comp_ptr is a pointer to a comparison function, e.g. std::less.
     cut.getter_ptr is a pointer to a member getter of SummaryEvent class. cut.value is a value
     the getter return is compared to. The lines in the loops expand to e.g.
     std::less<double>( evt->Get_track_energy(), 5 ). See below link for more info
     http://en.cppreference.com/w/cpp/utility/functional/function
   */
  
  // and cuts start with true, one failing cut fails all cuts
  Bool_t PassesAndCuts = true;

  for (auto cut: fAndCuts) {
    Bool_t PassesThisCut = cut.comp_ptr( cut.getter_ptr(*evt), cut.value );
    PassesAndCuts = PassesAndCuts && PassesThisCut;
  }

  // or cuts start with false, one true cut sets overall to true
  Bool_t PassesOrCuts = false;

  for (auto cut: fOrCuts) {
    Bool_t PassesThisCut = cut.comp_ptr( cut.getter_ptr(*evt), cut.value );
    PassesOrCuts = PassesOrCuts || PassesThisCut;
  }

  // if no or cuts are specified, set this to true, otherwise return and won't work
  if ( fOrCuts.size() == 0 ) PassesOrCuts = true;

  return (PassesAndCuts && PassesOrCuts);

}

//***************************************************************************************

/** 
    Function that copies variables from a `SummaryEvent` to member observables, depending on reco type.

    \param evt  Pointer to a summary event.

*/
void EventFilter::SetObservables(SummaryEvent *evt) {

  switch (fRecoType) {

  case mc_truth:
    fEnergy   = evt->Get_MC_energy();
    fBy       = evt->Get_MC_bjorkeny();
    fDir      = evt->Get_MC_dir();
    fPos      = evt->Get_MC_pos();
    break;
    
  case track:
    fEnergy   = evt->Get_track_energy();
    fBy       = evt->Get_track_bjorkeny();
    fDir      = evt->Get_track_dir();
    fPos      = evt->Get_track_pos();
    break;
    
  case shower:
    fEnergy   = evt->Get_shower_energy();
    fBy       = evt->Get_shower_bjorkeny();
    fDir      = evt->Get_shower_dir();
    fPos      = evt->Get_shower_pos();
    break;
   
  case customreco:
    
    if ( &GetCustomEnergy == NULL || &GetCustomDir == NULL || &GetCustomPos == NULL || GetCustomBy == NULL) {
      throw std::logic_error("ERROR! EventFilter::SetObservables() When using custom observables, pointers to the functions that return the observables need to be configured first by calling EventFilter::SetObsFuncPtrs(...)");
    }

    fEnergy = GetCustomEnergy(evt);
    fBy     = GetCustomBy(evt);
    fDir    = GetCustomDir(evt);
    fPos    = GetCustomPos(evt);

  }
  
}

//***************************************************************************************

/**
   This function has to be used to set the function pointers when reconstruction type "customreco" is used.

   Although it may seem complicated, the idea is actually rather simple. For example, consider that one wishes to use the direction of track reco, but the energy of a shower reco. This would be easy to implement without any pointers in the `switch` statement in `EventFilter::SetObservables`. However, typically one requires further logic, i.e. use shower reco energy only if it's quality level is at 1, otherwise use track reco. And then there is another thought that use shower energy estimate in one track score region and track energy estimate in another track score region. What I am trying to illustrate here, is that it will be unsustainable to add another switch element to `EventFilter::SetObservables` for each specific case.

   For this reason, an elegant alternative is provided. For such cases as illustrated above, the user can choose `customreco` at initialisation. Then, the user needs to define four functions in his/her application:
   \code{.cpp}
   Double_t MyCustomEnergy(SummaryEvent* evt) {...};
   TVector3 MyCustomDir(SummaryEvent* evt) {...};
   TVector3 MyCustomPos(SummaryEvent* evt) {...};
   Double_t MyCustomBY(SummaryEvent* evt) {...};
   \endcode
   In each function, the user has access to all of the data members of the `SummaryEvent` class (through the getter's) to make a specific selection which reconstruction variable is to be used in which circumstances. For example, if nothing specific is required for the direction reconstruction, the user can simply define:
   \code{.cpp}
   TVector3 MyCustomDir(SummaryEvent* evt) { return evt->Get_track_dir() };
   \endcode
   Having defined the four functions, the user needs to use this function to let the class know what to call. Example usage:
   \code{.cpp}
   EventFilter f1(EventFilter::customreco);
   f1.SetObsFuncPtrs( &MyCustomEnergy, &MyCustomDir, &MyCustomPos, &MyCustomBY);
   \endcode
   \param E   Address of the function that returns the custom energy
   \param ct  Address of the function that returns the custom direction vector
   \param pos Address of the function that returns the custom position vector
   \param by  Address of the function that returns the custom bjorken-y

 */
void EventFilter::SetObsFuncPtrs( Double_t (*E)(SummaryEvent*)  , TVector3 (*ct)(SummaryEvent*), 
				  TVector3 (*pos)(SummaryEvent*), Double_t (*by)(SummaryEvent*) ) {

  GetCustomEnergy = E;
  GetCustomDir = ct;
  GetCustomPos = pos;
  GetCustomBy  = by;

  std::cout << "NOTICE EventFilter::SetObsFuncPtrs() setting observable function pointers" << std::endl;

}
