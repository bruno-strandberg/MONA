#include "TH2.h"
#include "TRandom3.h"
#include <map>

/**
 * Namespace that holds functions and variables used in the GSGSampler.C macro.
 */
namespace GSGS {

  //forward declare struct evtid for function declarations
  struct evtid;

  //*****************************************************************
  // functions
  //*****************************************************************
  Bool_t   GetIntHists(TString flux_chain_file, Int_t flavor, Int_t is_cc);
  void     InitVars(Int_t flavor, Int_t is_cc);
  void     CleanUp();
  Double_t ReadGSGData(TString gsg_file_list, Int_t flavor, Int_t is_cc);
  void     CacheGSGdata(TString fname);
  Double_t ReadFromCache(TString fname);
  Bool_t   SampleEvents(TH2D *h_expected, TH2D *h_smeared,
			vector<evtid> **store, vector<evtid> **sample,
			Int_t low_sample_lim);
  void     StoreForWriting(Bool_t SampleOK, TH2D *smeared_nu, TH2D *smeared_nub, Int_t F, Int_t N);
  void     WriteToFiles(Int_t flavor, Int_t is_cc);
  TString  GetSummaryName(Int_t flavor, Int_t is_cc, Int_t emin, Int_t emax, Int_t runnr);

  //*****************************************************************
  // globally used variables in this script
  //*****************************************************************
  TH2D *fhInt_nu=NULL;   //!< histogram with expected number of interacted nu events
  TH2D *fhInt_nub=NULL;  //!< histogram with expected number of interacted nubar events
  TH2D *fhGSG_nu=NULL;   //!< histogram with all available MC nu events in GSG files
  TH2D *fhGSG_nub=NULL;  //!< histogram with all available MC nubar events in GSG files
  TRandom3 *fRand=NULL;  //!< random number generator
  Int_t    fEbins;       //!< number of energy bins in the event vectors
  Int_t    fCtbins;      //!< number of costheta bins in the event vectors
  Double_t fVcan = 0.;   //!< can size, set in ReadGSGData()
  Double_t fRhoSW = 0.;  //!< sea water density, set in ReadGSGData()

  //! 2D array of vectors, each vector holds summary nu eventid's of (E, costheta) bin
  vector<evtid> **fGSGEvts_nu;
  //! 2D array of vectors, each vector holds summary nub eventid's numbers of (E, costheta) bin
  vector<evtid> **fGSGEvts_nub;
  //! 2D array of vectors to hold a sub-sample of events in fGSGEvts_nu
  vector<evtid> **fSampleEvts_nu;
  //! 2D array of vectors to hold a sub-sample of events in fGSGEvts_nub
  vector<evtid> **fSampleEvts_nub;
  //! vector of vectors; each vector contains an event sample of one experiment
  vector< vector<evtid> > fExps;
  //! vector of vectors; each vector fExpHists[N] contains pointers to hists associated with fExps[N]
  vector< vector<TH2D*> > fExpHists;
  //! vector of strings; each string contains the flux file and sample index of an experiment
  vector< TString > fExpNames;

  //! map of flavor numbers and strings
  map < Int_t, TString > fFlavs  = { {0, "elec" },
				     {1, "muon" },
				     {2, "tau"  } };

  //! map of interaction numbers and strings
  map < Int_t, TString > fInts  = { {0, "NC" },
				    {1, "CC" } };

  /**
   *  Struct to store the Monte Carlo run number, event number and start of the energy range.
   */
  struct evtid {
    Int_t run_nr; //!< gSeaGen run number
    Int_t evt_nr; //!< gSeaGen event number
    Int_t e_min;  //!< Minimum energy of the neutrinos in gSeaGen file (1 or 3)

    //! Default constructor
    evtid(): run_nr(0), evt_nr(0), e_min(0) {}
    
    /**
     * Constructor
     * \param _run_nr gSeaGen file number
     * \param _evt_nr gSeaGen event number in file
     * \param _e_min  Minimum neutrino energy in gSeaGen file
     */
    evtid(Int_t _run_nr, Int_t _evt_nr, Int_t _e_min) {
      run_nr = _run_nr;
      evt_nr = _evt_nr;
      e_min  = _e_min;
    }

    /**
     * Operator to sort events by increasing file number, event id and energy by using std::sort
     * \param   rhs instance of evtid that this object is compared to
     * \return  true if second run nr is larger; if equal run numbers true if second e_min larger;
     *          if equal run nr and equal e_min true if second event number larger.
     */
    bool operator < (const evtid& rhs) const {
 
      if ( this->run_nr < rhs.run_nr ) { 
	return true;         
      }

      else if ( this->run_nr == rhs.run_nr ) {

	if ( this->e_min < rhs.e_min) {
	  return true;
	}
	else if (this->e_min == rhs.e_min) {
	  return (this->evt_nr < rhs.evt_nr); 
	}
	else {
	  return false;
	}

      }

      else {
	return false; 
      }

    }

    /** 
     * Operator to find identical events by using std::find.
     * \param rhs instance of evtid that this object is compared to
     * \return true if run_nr and evt_nr and e_min of this object and other object match.
     */
    bool operator == (const evtid& rhs) const {
      return ( ( this->run_nr == rhs.run_nr ) && 
	       ( this->evt_nr == rhs.evt_nr ) && 
	       ( this->e_min  == rhs.e_min  ) );
    }

    /** 
     * Not equal operator.
     * \param rhs instance of evtid that this object is compared to
     * \return true if run_nr or evt_nr or e_min of this object and other object mismatch.
     */
    bool operator != (const evtid &rhs) const {
      return !(*this == rhs);
    }

    /**
     * Stream operator for cout.
     */
    friend std::ostream &operator << ( std::ostream &output, const evtid &EID ) { 
      output << EID.run_nr << ' ' << EID.evt_nr << ' ' << EID.e_min;  
      return output;  
    }

  };

  /** 
   * Function to find the first and last index of the elements that correspond to the same gSeaGen
   * file in a sorted vector of evtid's.
   * .
   *
   * \param evts    Sorted vector of evtid's.
   * \param _run_nr gSeaGen run number.
   * \param _e_min  gSeaGen energy range start.
   * \return Returns a pair of indices; evts[pair.first] is the location of the first
   *         evtid in vector evts with run number _run_nr and energy start _e_min; evts[pair.second]
   *         is the location of the last evtid with run number _run_nr and energy start _e_min.
   *  
   */
  std::pair<Int_t, Int_t> get_start_stop(vector<evtid> evts, Int_t _run_nr, Int_t _e_min) {

    Int_t start = evts.size();
    Int_t stop  = evts.size();

    for (Int_t i = 0; i < (Int_t)evts.size(); i++) {

      if (evts[i].run_nr == _run_nr && evts[i].e_min == _e_min) {

	start = i;
	stop  = start;

	while ( evts[stop].run_nr == _run_nr && evts[stop].e_min == _e_min ) stop++;

	break;

      }

    } //end for

    return std::make_pair(start, stop);

  }

}
