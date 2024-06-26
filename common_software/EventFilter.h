#ifndef EventFilter_h
#define EventFilter_h
#include "SummaryEvent.h"

#include<vector>
#include<map>
#include <functional>

/**
   Structure for use in `EventFilter` that stores info necessary to define a cut on a `SummaryEvent`.
*/
struct cutobj {

  std::function< Double_t(SummaryEvent&) > getter_ptr; //!< pointer to a SummaryEvent getter func
  Double_t value;                                      //!< value the getter return is compared to
  std::function< bool(double, double) > comp_ptr;      //!< pointer to comparison function

  /** Default constructor */
cutobj(): getter_ptr(0), value(0), comp_ptr(0) {};
    
  /**
     Constructor
     \param _getter_ptr Pointer to a SummaryEvent getter function
     \param _value      Value the getter function return is compared to
     \param _comp_ptr   Pointer to an std comparison func
  */
  cutobj(std::function<Double_t(SummaryEvent&)> _getter_ptr, Double_t _value,
	 std::function<bool(double, double)> _comp_ptr) {
    getter_ptr = _getter_ptr;
    value      = _value;
    comp_ptr   = _comp_ptr;
  }

};

/**
   A class to filter events.

   The class implements an event filter for the event class `SummaryEvent`. The function `AddCut`
   can be used to define event selections, the function `PassesCuts(SummaryEvent *evt)` returns 
   true/false, depending on whether the argument event passed the selection or not.

   Example usage:
   \code{.cpp}
   EventFilter f;
   f.AddCut( &SummaryEvent::Get_track_dir_z , std::less<double>()         ,   0, true );
   f.AddCut( &SummaryEvent::Get_track_energy, std::greater_equal<double>(),   5, true );
   f.AddCut( &SummaryEvent::Get_shower_ql1  , std::greater_equal<double>(), 0.5, false );
   f.AddCut( &SummaryEvent::Get_track_ql1   , std::equal_to<double>()     ,   1, false );
   \endcode
   In this case, the function `PassesCuts(SummaryEvent *evt)` will return
   \code{.cpp}
   ( ( evt->GetTrack_dir_z() < 0 ) && ( evt->Get_track_energy() >= 5 ) ) && 
   ( ( evt->Get_shower_ql1() >= 0.5 ) || ( evt->Get_track_ql1() == 0.5 ) )
   \endcode

   As evident from the example, the function `AddCut` takes as arguments a pointer to a getter 
   function of the `SummaryEvent` class, a value that the return of the getter function
   is compared to, a pointer to a comparison operator and a flag to indicate whether this cut is an 
   'and' cut or an 'or' cut. As long as only double's are used in `SummaryEvent`, the standard library
   functions `std::less<double>()` , ... etc are the comparison operators to pass.

   Additionally, the class has variables for observables energy, position, direction and bjorken-y, and
   a variable to determine which reconstruction type (e.g. track, shower, mc_truth) is used to set the
   observables. The observables are typically used in inheriting classes (`EventSelection` and 
   `DetResponse`) in filling data structures.

   There is a slightly more advanced option to use a custom reconstruction type. This is very useful in cases where, say, the user wishes to use shower energy estimate with track direction estimate. Consult the documentation of the function `EventFilter::SetObsFuncPtrs` for details how to achieve this.

   Examples and test for `EventFilter` and inheriting classes `DetResponse` and `EventSelection` can be found at `NMH/tests/common_software`.

   Future developments: if SummaryEvent is changed to consist also of other variable types than
   doubles, then this class needs to be templated.
 */
class EventFilter {

 public:

  /// enum to determine which reco type is used in filling hists in inheriting classes
  enum reco {
    mc_truth,  //!< use monte-carlo truth for energy, direction, bjorken-y and position
    track,     //!< use track-reco for energy, direction, bjorken-y and position
    shower,    //!< use shower-reco for energy, direction, bjorken-y and position
    customreco //!< use user-defined energy, dir, bjorken-y and pos; see `EventFilter::SetObsFuncPtrs` doc.
  };
  
  EventFilter(reco reco_type = EventFilter::mc_truth);
  ~EventFilter();
  
  Bool_t PassesCuts(SummaryEvent *evt);
  void   SetObservables(SummaryEvent *evt);
  void   AddCut( std::function<Double_t(SummaryEvent&)> getter_ptr,
		 std::function<bool(double, double)> comp_ptr, Double_t value, Bool_t AndCut);
  void   AddCut( TString getter_func_name, TString comp_op, Double_t value, Bool_t AndCut);
  
  void   SetObsFuncPtrs( Double_t (*E)(SummaryEvent*)  , TVector3 (*ct)(SummaryEvent*), 
			 TVector3 (*pos)(SummaryEvent*), Double_t (*by)(SummaryEvent*) );

  reco     fRecoType;    //!< variable to determine which reco type is used
  Double_t fEnergy;      //!< energy observable for inheriting classes, depending on reco type
  TVector3 fPos;         //!< position observable for inheriting classes, depending on reco type
  TVector3 fDir;         //!< direction observable for inheriting classes, depending on reco type
  Double_t fBy;          //!< bjorken-y observable for inheriting classes, depending on reco type
  
 private:
  
  std::vector<cutobj> fAndCuts; //!< vector with 'and' cuts
  std::vector<cutobj> fOrCuts;  //!< vector with 'or' cuts

  Double_t (*GetCustomEnergy)(SummaryEvent*);
  TVector3 (*GetCustomDir)(SummaryEvent*);
  TVector3 (*GetCustomPos)(SummaryEvent*);
  Double_t (*GetCustomBy)(SummaryEvent*);

  /// map of comparison characters and comparison operators, required to be able to use MONA with pyroot
  std::map< TString, std::function<bool(double, double)> > fOpMap = {
    { ">" , std::greater<double>()       },
    { ">=", std::greater_equal<double>() },
    { "==", std::equal_to<double>()      },
    { "<=", std::less_equal<double>()    },
    { "<" , std::less<double>()          },
  };
    
};

#endif
