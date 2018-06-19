#ifndef EventFilter_h
#define EventFilter_h
#include "SummaryEvent.h"
#include<vector>

/**
   A class to filter events.

   The class has two functions to define an event filter for the event class SummaryEvent.
   The functions `AddAndCut` and `AddOrCut` can be used to define event selections, the function
   `PassesCuts(SummaryEvent *evt)` returns true/false, depending on whether the argument event
   passed the selection or not.

   Example usage:
   ```
   EventFilter f;
   f.AddAndCut( &SummaryEvent::Get_track_dir_z ,   0, std::less<double>() );
   f.AddAndCut( &SummaryEvent::Get_track_energy,   5, std::greater_equal<double>() );
   f.AddOrCut(  &SummaryEvent::Get_shower_ql1  , 0.5, std::greater_equal<double>() );
   f.AddOrCut(  &SummaryEvent::Get_track_ql1   ,   1, std::equal<double>() );
   ```
   In this case, the function `PassesCuts(SummaryEvent *evt)` will return
   ```
   ( ( evt->GetTrack_dir_z() < 0 ) && ( evt->Get_track_energy() >= 5 ) ) && 
   ( ( evt->Get_shower_ql1() >= 0.5 ) || ( evt->Get_track_ql1() == 0.5 ) )
   ```

   As evident from the example, the functions `AddAndCut`, `AddOrCut` take as arguments a pointer
   to a getter function of the SummaryEvent class, a value that the return of the getter function
   is compared to, and a pointer to a comparison operator. As long as only double's are used in
   SummaryEvent, the standard library functions std::less<double>() , ... etc are the comparison
   operators to pass.

   Future developments: if SummaryEvent is changed to consist also of other variable types than
   doubles, then this class needs to be templated.
 */
class EventFilter {

 public:
  EventFilter() {};
  ~EventFilter() {};
  
  Bool_t PassesCuts(SummaryEvent *evt);
  void   AddAndCut( std::function<Double_t(SummaryEvent&)> getter_ptr,
		    Double_t value,
		    std::function<bool(double, double)> comp_ptr);
  void   AddOrCut( std::function<Double_t(SummaryEvent&)> getter_ptr,
		   Double_t value,
		   std::function<bool(double, double)> comp_ptr);

 private:

  /**
     Structure for internal use that stores info necessary to define a cut on a SummaryEvent object
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

  std::vector<cutobj> fAndCuts; //!< vector with 'and' cuts
  std::vector<cutobj> fOrCuts;  //!< vector with 'or' cuts

};

#endif
