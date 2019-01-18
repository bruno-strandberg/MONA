#include "EventSelection.h"
#include "SummaryParser.h"
#include "TVector3.h"
#include "TMath.h"

/** 
    Namespace for functions of custom EventFilter. 
*/
namespace CUSTOMEF {

  /** Function to return hybrid energy, used instead of track or shower energy
      \param evt  Pointer to a summary event instance
      \return     reconstructed energy
   */
  Double_t HybridEnergy(SummaryEvent *evt) {
    // A misreconstructed shower has a negative energy.
    if ((evt->Get_shower_energy()) > 0 and (evt->Get_shower_energy() <= 100)) {
      return evt->Get_shower_energy();
    }
    else {
      return evt->Get_track_energy();
    }
  }

  // for direction, position and bjorken-y I just use the track standards
  TVector3 TrackDir(SummaryEvent* evt) { return evt->Get_track_dir();      }
  TVector3 TrackPos(SummaryEvent* evt) { return evt->Get_track_pos();      }
  Double_t TrackBY (SummaryEvent* evt) { return evt->Get_track_bjorkeny(); }

};


