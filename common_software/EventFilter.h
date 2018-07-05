#ifndef EventFilter_h
#define EventFilter_h
#include "SummaryEvent.h"
#include<vector>

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
   ```
   EventFilter f;
   f.AddCut( &SummaryEvent::Get_track_dir_z , std::less<double>()         ,   0, true );
   f.AddCut( &SummaryEvent::Get_track_energy, std::greater_equal<double>(),   5, true );
   f.AddCut( &SummaryEvent::Get_shower_ql1  , std::greater_equal<double>(), 0.5, false );
   f.AddCut( &SummaryEvent::Get_track_ql1   , std::equal_to<double>()     ,   1, false );
   ```
   In this case, the function `PassesCuts(SummaryEvent *evt)` will return
   ```
   ( ( evt->GetTrack_dir_z() < 0 ) && ( evt->Get_track_energy() >= 5 ) ) && 
   ( ( evt->Get_shower_ql1() >= 0.5 ) || ( evt->Get_track_ql1() == 0.5 ) )
   ```

   As evident from the example, the function `AddCut` takes as arguments a pointer to a getter 
   function of the `SummaryEvent` class, a value that the return of the getter function
   is compared to, a pointer to a comparison operator and a flag to indicate whether this cut is an 
   'and' cut or an 'or' cut. As long as only double's are used in `SummaryEvent`, the standard library
   functions `std::less<double>()` , ... etc are the comparison operators to pass.

   Future developments: if SummaryEvent is changed to consist also of other variable types than
   doubles, then this class needs to be templated.
 */
class EventFilter {

 public:
  EventFilter() {};
  ~EventFilter() {};
  
  Bool_t PassesCuts(SummaryEvent *evt);
  void   AddCut( std::function<Double_t(SummaryEvent&)> getter_ptr,
		 std::function<bool(double, double)> comp_ptr, Double_t value, Bool_t AndCut);

 private:

  std::vector<cutobj> fAndCuts; //!< vector with 'and' cuts
  std::vector<cutobj> fOrCuts;  //!< vector with 'or' cuts

};

#endif
