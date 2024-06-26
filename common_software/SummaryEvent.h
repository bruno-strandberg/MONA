#include "TObject.h"
#include "TVector3.h"
#include "TRandom3.h"
#include <functional>
#include <map>

#ifndef SummaryEvent_h
#define SummaryEvent_h

/**
   Data format for NMH analysis.

   It is obvious that for a code to be maintainable, one requires a certain format for the data. At the time when this code was started, there was no clear agreement on what is the data format that is input to the NMH sensitivity analyses. For this reason this event class was defined, the data from ECAP Random Decision Forest PID is converted to this format by scripts in `NMH/data_sorting/`.

   SummaryEvent consists of double's only. The main reason for this is to fascilitate simple data input to RooFit (where NMH fitting is performed), which works with flat trees. Additionally, such a tree can be easilty inspected/accessed without access to the SummaryEvent class. The variables are all private to dis-allow access without the use of the interface (setters/getters). The existing function names and variable names should not be changed, but new ones can be added.

   The NMH sensitivy is itself a very high-level analysis. For this reason there are placeholders only for one track reco and one shower reco. It is envisaged that the selection of the 'best' fitted track is done before this analysis and hence this data format does not include several tracks per event as e.g. aanet.

   Similarly, currently there are placeholders for three levels of quality cuts, more can be added. Quality levels (i.e. event selection cuts) should be defined prior to this analysis or when data is written to this format. Example: one may wish to have all reconstructed track events (ql0), reconstructed track events with likelihood higher than some value (ql1), reconstructed tracks with likelihood higher than some value and a certain number of hits (ql2), etc. 

   If a better reco (some DNN) or PID becomes available, this info can either be stored in the existing variables or, if necessary, new variables can be added.
 */
class SummaryEvent : public TObject {

 public:
  SummaryEvent();
  ~SummaryEvent();

  void SetSeed(UInt_t seed);
  
  //public interface to the data structure

  // setters
  void Set_MC_runID(Double_t mc_runid)                { fMC_runID = mc_runid; }       
  void Set_MC_evtID(Double_t mc_evtid)                { fMC_evtID = mc_evtid; }       
  void Set_MC_w2(Double_t mc_w2)                      { fMC_w2 = mc_w2; }          
  void Set_MC_w1y(Double_t mc_w1y)                    { fMC_w1y = mc_w1y; }
  void Set_MC_w2denom(Double_t mc_w2denom)            { fMC_w2denom = mc_w2denom; }       
  void Set_MC_erange_start(Double_t mc_Estart)        { fMC_erange_start = mc_Estart; }
  void Set_MC_erange_stop(Double_t mc_Estop)          { fMC_erange_stop = mc_Estop; }
  void Set_MC_is_CC(Double_t mc_iscc)                 { fMC_is_CC = mc_iscc; }  
  void Set_MC_is_neutrino(Double_t mc_isnu)           { fMC_is_neutrino = mc_isnu; }
  void Set_MC_type(Double_t mc_type)                  { fMC_type = mc_type; }
  void Set_MC_ichan(Double_t mc_ichan)                { fMC_ichan = mc_ichan; }
  void Set_MC_energy(Double_t mc_e)                   { fMC_energy = mc_e; }
  void Set_MC_bjorkeny(Double_t mc_by)                { fMC_bjorkeny = mc_by; }
  void Set_MC_dir(Double_t x, Double_t y, Double_t z) { 
    fMC_dir_x = x; 
    fMC_dir_y = y; 
    fMC_dir_z = z;
  }
  void Set_MC_pos(Double_t x, Double_t y, Double_t z) { 
    fMC_pos_x = x; 
    fMC_pos_y = y;
    fMC_pos_z = z; 
  }

  void Set_track_energy(Double_t t_e)                 { fTrack_energy = t_e; }
  void Set_track_bjorkeny(Double_t t_by)              { fTrack_bjorkeny = t_by; }
  void Set_track_ql0(Double_t t_ql0)                  { fTrack_ql0 = t_ql0; }     
  void Set_track_ql1(Double_t t_ql1)                  { fTrack_ql1 = t_ql1; }
  void Set_track_ql2(Double_t t_ql2)                  { fTrack_ql2 = t_ql2; }
  void Set_track_ql3(Double_t t_ql3)                  { fTrack_ql3 = t_ql3; }
  void Set_track_dir(Double_t x, Double_t y, Double_t z) { 
    fTrack_dir_x = x; 
    fTrack_dir_y = y; 
    fTrack_dir_z = z;
  }
  void Set_track_pos(Double_t x, Double_t y, Double_t z) { 
    fTrack_pos_x = x; 
    fTrack_pos_y = y;
    fTrack_pos_z = z;
  }

  void Set_shower_energy(Double_t s_e)                { fShower_energy = s_e; }
  void Set_shower_bjorkeny(Double_t s_by)             { fShower_bjorkeny = s_by; }
  void Set_shower_ql0(Double_t s_ql0)                 { fShower_ql0 = s_ql0; }     
  void Set_shower_ql1(Double_t s_ql1)                 { fShower_ql1 = s_ql1; }
  void Set_shower_ql2(Double_t s_ql2)                 { fShower_ql2 = s_ql2; }
  void Set_shower_ql3(Double_t s_ql3)                 { fShower_ql3 = s_ql3; }
  void Set_shower_dir(Double_t x, Double_t y, Double_t z) { 
    fShower_dir_x = x; 
    fShower_dir_y = y; 
    fShower_dir_z = z;
  }
  void Set_shower_pos(Double_t x, Double_t y, Double_t z) { 
    fShower_pos_x = x; 
    fShower_pos_y = y; 
    fShower_pos_z = z;
  }

  void Set_RDF_muon_score(Double_t rms)               { fRDF_muon_score = rms; } 
  void Set_RDF_track_score(Double_t rts)              { fRDF_track_score = rts; }
  void Set_RDF_noise_score(Double_t rns)              { fRDF_noise_score = rns; }

  // getters
  Double_t Get_MC_runID()        { return fMC_runID; }       
  Double_t Get_MC_evtID()        { return fMC_evtID; }       
  Double_t Get_MC_w2()           { return fMC_w2; }          
  Double_t Get_MC_w1y()          { return fMC_w1y; }
  Double_t Get_MC_w2denom()      { return fMC_w2denom; }       
  Double_t Get_MC_erange_start() { return fMC_erange_start; }
  Double_t Get_MC_erange_stop()  { return fMC_erange_stop;  }
  Double_t Get_MC_is_CC()        { return fMC_is_CC; }  
  Double_t Get_MC_is_neutrino()  { return fMC_is_neutrino; }
  Double_t Get_MC_type()         { return fMC_type; }
  Double_t Get_MC_ichan()        { return fMC_ichan; }
  Double_t Get_MC_dir_x()        { return fMC_dir_x; }
  Double_t Get_MC_dir_y()        { return fMC_dir_y; }
  Double_t Get_MC_dir_z()        { return fMC_dir_z; }
  Double_t Get_MC_pos_x()        { return fMC_pos_x; }
  Double_t Get_MC_pos_y()        { return fMC_pos_y; }
  Double_t Get_MC_pos_z()        { return fMC_pos_z; }
  Double_t Get_MC_energy()       { return fMC_energy; }
  Double_t Get_MC_bjorkeny()     { return fMC_bjorkeny; }
  TVector3 Get_MC_dir()          { return TVector3(fMC_dir_x, fMC_dir_y, fMC_dir_z); }
  TVector3 Get_MC_pos()          { return TVector3(fMC_pos_x, fMC_pos_y, fMC_pos_z); }

  Double_t Get_track_dir_x()     { return fTrack_dir_x; }
  Double_t Get_track_dir_y()     { return fTrack_dir_y; }
  Double_t Get_track_dir_z()     { return fTrack_dir_z; }
  Double_t Get_track_pos_x()     { return fTrack_pos_x; }
  Double_t Get_track_pos_y()     { return fTrack_pos_y; }
  Double_t Get_track_pos_z()     { return fTrack_pos_z; }
  Double_t Get_track_energy()    { return fTrack_energy; }
  Double_t Get_track_bjorkeny()  { return fTrack_bjorkeny; }
  Double_t Get_track_ql0()       { return fTrack_ql0; }     
  Double_t Get_track_ql1()       { return fTrack_ql1; }
  Double_t Get_track_ql2()       { return fTrack_ql2; }
  Double_t Get_track_ql3()       { return fTrack_ql3; }
  TVector3 Get_track_dir()       { return TVector3(fTrack_dir_x, fTrack_dir_y, fTrack_dir_z); }
  TVector3 Get_track_pos()       { return TVector3(fTrack_pos_x, fTrack_pos_y, fTrack_pos_z); }

  Double_t Get_shower_dir_x()    { return fShower_dir_x; }
  Double_t Get_shower_dir_y()    { return fShower_dir_y; }
  Double_t Get_shower_dir_z()    { return fShower_dir_z; }
  Double_t Get_shower_pos_x()    { return fShower_pos_x; }
  Double_t Get_shower_pos_y()    { return fShower_pos_y; }
  Double_t Get_shower_pos_z()    { return fShower_pos_z; }
  Double_t Get_shower_energy()   { return fShower_energy; }
  Double_t Get_shower_bjorkeny() { return fShower_bjorkeny; }
  Double_t Get_shower_ql0()      { return fShower_ql0; }     
  Double_t Get_shower_ql1()      { return fShower_ql1; }
  Double_t Get_shower_ql2()      { return fShower_ql2; }
  Double_t Get_shower_ql3()      { return fShower_ql3; }
  TVector3 Get_shower_dir()      { return TVector3(fShower_dir_x, fShower_dir_y, fShower_dir_z); }
  TVector3 Get_shower_pos()      { return TVector3(fShower_pos_x, fShower_pos_y, fShower_pos_z); }

  Double_t Get_RDF_muon_score()  { return fRDF_muon_score; } 
  Double_t Get_RDF_track_score() { return fRDF_track_score; }
  Double_t Get_RDF_noise_score() { return fRDF_noise_score; }
  
  // other functions
  void FillPseudoData(Bool_t logE=kFALSE);

 private:

  // variables of SummaryEvent; this is designed as a flat tree of Double_t's for
  // easier interfacing with RooFit (does not like class structures as branches)

  Double_t        fMC_runID;        //!< gSeaGen file run number
  Double_t        fMC_evtID;        //!< Event ID in gSeaGen run
  Double_t        fMC_w2;           //!< Weight 2 from gSeaGen
  Double_t        fMC_w1y;          //!< multiply atm. muon/noise/nu by this to get evts in year. Calculated as part of ECAP PID pipeline, valid if all events from production are used
  Double_t        fMC_w2denom;      //!< number of generated events in gSeaGen run; for atm. muons and noise the livetime of the run in seconds
  Double_t        fMC_erange_start; //!< Start of gSeaGen E range, e.g. for 3--100 GeV this is 3
  Double_t        fMC_erange_stop;  //!< Stop of gSeaGen E range, e.g. for 3--100 GeV this is 100
  Double_t        fMC_is_CC;        //!< 1 - CC, 0 - NC
  Double_t        fMC_is_neutrino;  //!< 1 - atm neutrion, 0 - atm muon
  Double_t        fMC_type;         //!< PDG code of the primary lepton
  Double_t        fMC_ichan;        //!< neutrino interaction channel from gSeaGen
  Double_t        fMC_dir_x;
  Double_t        fMC_dir_y;
  Double_t        fMC_dir_z;
  Double_t        fMC_pos_x;
  Double_t        fMC_pos_y;
  Double_t        fMC_pos_z;
  Double_t        fMC_energy;
  Double_t        fMC_bjorkeny;

  Double_t        fTrack_dir_x;
  Double_t        fTrack_dir_y;
  Double_t        fTrack_dir_z;
  Double_t        fTrack_pos_x;
  Double_t        fTrack_pos_y;
  Double_t        fTrack_pos_z;
  Double_t        fTrack_energy;    //!< Reconstructed neutrino energy
  Double_t        fTrack_bjorkeny;  //!< Reconstructed Bjorken Y (currently only placeholder)
  Double_t        fTrack_ql0;       //!< Quality level 0 - lowest quality event (few quality checks), typically just that reco worked
  Double_t        fTrack_ql1;       //!< Quality level 1 - more quality checks, user can select what to store here, see `MONA/README.md` section `Quality levels`
  Double_t        fTrack_ql2;       //!< Quality level 2 - even more quality checks, user can select what to store here, see `MONA/README.md` section `Quality levels`
  Double_t        fTrack_ql3;       //!< Quality level 3 - even more quality checks, user can select what to store here, see `MONA/README.md` section `Quality levels`

  Double_t        fShower_dir_x;
  Double_t        fShower_dir_y;
  Double_t        fShower_dir_z;
  Double_t        fShower_pos_x;
  Double_t        fShower_pos_y;
  Double_t        fShower_pos_z;
  Double_t        fShower_energy;   //!< Reconstructed neutrino energy
  Double_t        fShower_bjorkeny; //!< Reconstructed Bjorken Y
  Double_t        fShower_ql0;      //!< Similar to gandalf quality levels, see `MONA/README.md` section `Quality levels`
  Double_t        fShower_ql1;
  Double_t        fShower_ql2;
  Double_t        fShower_ql3;

  Double_t        fRDF_muon_score;  //!< PID score for this being an atm muon from random dec. forest
  Double_t        fRDF_track_score; //!< PID score for this being a track event from random dec. forest
  Double_t        fRDF_noise_score; //!< PID score for this being a noise event from random dec. forest

  // class definition and version. The index is added to the output. If you e.g. add some more
  // variables, increment the intex by 1. Then at read-in you can check the index and e.g. ignore
  // some variables in older versions
  ClassDef(SummaryEvent, 2)
};

#endif

#ifndef SUMMARYUTILS
#define SUMMARYUTILS

/** This namespace provides some structures and functions for miscellaneous operations involving `SummaryEvent` data. */
namespace SUMMARYUTILS {

  TRandom3 fRand(0); //!< random generator for `SummaryEvent::FillPseudoData`
  std::vector<Double_t> fFlavs = {12.,14.,16.,-12.,-14.,-16.}; //!< PDG codes for `SummaryEvent::FillPseudoData`

  /// map of function names and corresponding getter functions from `SummaryEvent` class, this is required to be able to use MONA with PyROOT. If more variables & getters are added to `SummaryEvent`, the getter functions should be added to this map to make the new functionality also available in PyROOT. The map is used in the `EventFilter::AddCut(TString, TString, Double_t, Bool_t)` function.
  
  std::map< TString, std::function<Double_t(SummaryEvent&)> > fFuncMap = {
    {"SummaryEvent.Get_MC_runID"       , &SummaryEvent::Get_MC_runID       },
    {"SummaryEvent.Get_MC_evtID"       , &SummaryEvent::Get_MC_evtID       },
    {"SummaryEvent.Get_MC_w2"          , &SummaryEvent::Get_MC_w2          },
    {"SummaryEvent.Get_MC_w1y"         , &SummaryEvent::Get_MC_w1y         },
    {"SummaryEvent.Get_MC_w2denom"     , &SummaryEvent::Get_MC_w2denom     },  
    {"SummaryEvent.Get_MC_erange_start", &SummaryEvent::Get_MC_erange_start},
    {"SummaryEvent.Get_MC_erange_stop" , &SummaryEvent::Get_MC_erange_stop },
    {"SummaryEvent.Get_MC_is_CC"       , &SummaryEvent::Get_MC_is_CC       },
    {"SummaryEvent.Get_MC_is_neutrino" , &SummaryEvent::Get_MC_is_neutrino },
    {"SummaryEvent.Get_MC_type"        , &SummaryEvent::Get_MC_type        },
    {"SummaryEvent.Get_MC_ichan"       , &SummaryEvent::Get_MC_ichan       },
    {"SummaryEvent.Get_MC_dir_x"       , &SummaryEvent::Get_MC_dir_x       },
    {"SummaryEvent.Get_MC_dir_y"       , &SummaryEvent::Get_MC_dir_y       },
    {"SummaryEvent.Get_MC_dir_z"       , &SummaryEvent::Get_MC_dir_z       },
    {"SummaryEvent.Get_MC_pos_x"       , &SummaryEvent::Get_MC_pos_x       },
    {"SummaryEvent.Get_MC_pos_y"       , &SummaryEvent::Get_MC_pos_y       },
    {"SummaryEvent.Get_MC_pos_z"       , &SummaryEvent::Get_MC_pos_z       },
    {"SummaryEvent.Get_MC_energy"      , &SummaryEvent::Get_MC_energy      },
    {"SummaryEvent.Get_MC_bjorkeny"    , &SummaryEvent::Get_MC_bjorkeny    },
    {"SummaryEvent.Get_track_dir_x"    , &SummaryEvent::Get_track_dir_x    },
    {"SummaryEvent.Get_track_dir_y"    , &SummaryEvent::Get_track_dir_y    },
    {"SummaryEvent.Get_track_dir_z"    , &SummaryEvent::Get_track_dir_z    },
    {"SummaryEvent.Get_track_pos_x"    , &SummaryEvent::Get_track_pos_x    },
    {"SummaryEvent.Get_track_pos_y"    , &SummaryEvent::Get_track_pos_y    },
    {"SummaryEvent.Get_track_pos_z"    , &SummaryEvent::Get_track_pos_z    },
    {"SummaryEvent.Get_track_energy"   , &SummaryEvent::Get_track_energy   },
    {"SummaryEvent.Get_track_bjorkeny" , &SummaryEvent::Get_track_bjorkeny },
    {"SummaryEvent.Get_track_ql0"      , &SummaryEvent::Get_track_ql0      },
    {"SummaryEvent.Get_track_ql1"      , &SummaryEvent::Get_track_ql1      },
    {"SummaryEvent.Get_track_ql2"      , &SummaryEvent::Get_track_ql2      },
    {"SummaryEvent.Get_track_ql3"      , &SummaryEvent::Get_track_ql3      },
    {"SummaryEvent.Get_shower_dir_x"   , &SummaryEvent::Get_shower_dir_x   },
    {"SummaryEvent.Get_shower_dir_y"   , &SummaryEvent::Get_shower_dir_y   },
    {"SummaryEvent.Get_shower_dir_z"   , &SummaryEvent::Get_shower_dir_z   },
    {"SummaryEvent.Get_shower_pos_x"   , &SummaryEvent::Get_shower_pos_x   },
    {"SummaryEvent.Get_shower_pos_y"   , &SummaryEvent::Get_shower_pos_y   },
    {"SummaryEvent.Get_shower_pos_z"   , &SummaryEvent::Get_shower_pos_z   },
    {"SummaryEvent.Get_shower_energy"  , &SummaryEvent::Get_shower_energy  },
    {"SummaryEvent.Get_shower_bjorkeny", &SummaryEvent::Get_shower_bjorkeny},
    {"SummaryEvent.Get_shower_ql0"     , &SummaryEvent::Get_shower_ql0     },
    {"SummaryEvent.Get_shower_ql1"     , &SummaryEvent::Get_shower_ql1     },
    {"SummaryEvent.Get_shower_ql2"     , &SummaryEvent::Get_shower_ql2     },
    {"SummaryEvent.Get_shower_ql3"     , &SummaryEvent::Get_shower_ql3     },
    {"SummaryEvent.Get_RDF_muon_score" , &SummaryEvent::Get_RDF_muon_score },
    {"SummaryEvent.Get_RDF_track_score", &SummaryEvent::Get_RDF_track_score},
    {"SummaryEvent.Get_RDF_noise_score", &SummaryEvent::Get_RDF_noise_score},
  };

};

#endif
