#include "EventFilter.h"

//***************************************************************************************

/**
   Function to add an 'and' cut to this event filter.

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
