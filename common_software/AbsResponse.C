#include "AbsResponse.h"
#include "NMHUtils.h"

/** Constructor with same binning settings for true and reco space.
    
    \param reco_type         Type of reco variables used to fill the response, e.g. DetResponse::track; See `EventFilter`
    \param resp_name         Name of the response
    \param ebins             Number of energy bins
    \param emin              Energy minimum, should match the minimum energy of simulated (gSeaGen) events
    \param emax              Energy maximum, should match the maximum energy of simulated (gSeaGen) events
    \param ctbins            Number of cos-theta bins
    \param ctmin             cos-theta minimum, should match the minimum cos-theta of simulated (gSeaGen) events
    \param ctmax             cos-theta maximum, should match the maximum cos-theta of simulated (gSeaGen) events
    \param bybins            Number of bjorken-y bins
    \param bymin             bjorken-y minimum, should match the minimum bjorken-y of simulated (gSeaGen) events
    \param bymax             bjorken-y maximum, should match the maximum bjorken-y of simulated (gSeaGen) events

 */
AbsResponse::AbsResponse(reco reco_type, TString resp_name,
			 Int_t ebins , Double_t emin , Double_t emax ,
			 Int_t ctbins, Double_t ctmin, Double_t ctmax,
			 Int_t bybins, Double_t bymin, Double_t bymax ) :
  AbsResponse(reco_type, resp_name, ebins, emin, emax, ctbins, ctmin, ctmax, bybins, bymin, bymax, ebins, emin, emax, ctbins, ctmin, ctmax, bybins, bymin, bymax) {}

//=========================================================================================================

/** Constructor with same binning settings for true and reco space.
    
    \param reco_type         Type of reco variables used to fill the response, e.g. DetResponse::track; See `EventFilter`
    \param resp_name         Name of the response
    \param h_bins            Pointer to a TH3 with desired bins, energy on the x-axis, cos-theta on the y-axis, bjorken-y on the z-axis
 */
AbsResponse::AbsResponse(reco reco_type, TString resp_name, TH3* h_bins) :
  AbsResponse(reco_type, resp_name, h_bins, h_bins) {}

//=========================================================================================================

/** Constructor with different binning settings for true and reco space.
    
    \param reco_type         Type of reco variables used to fill the response, e.g. DetResponse::track; See `EventFilter`
    \param resp_name         Name of the response

    \param t_ebins           Number of energy bins, true space
    \param t_emin            Energy minimum, should match the minimum energy of simulated (gSeaGen) events, true space
    \param t_emax            Energy maximum, should match the maximum energy of simulated (gSeaGen) events, true space
    \param t_ctbins          Number of cos-theta bins, true space
    \param t_ctmin           cos-theta minimum, should match the minimum cos-theta of simulated (gSeaGen) events, true space
    \param t_ctmax           cos-theta maximum, should match the maximum cos-theta of simulated (gSeaGen) events, true space
    \param t_bybins          Number of bjorken-y bins, true space
    \param t_bymin           bjorken-y minimum, should match the minimum bjorken-y of simulated (gSeaGen) events, true space
    \param t_bymax           bjorken-y maximum, should match the maximum bjorken-y of simulated (gSeaGen) events, true space

    \param r_ebins           Number of energy bins, reco space
    \param r_emin            Energy minimum, reco space
    \param r_emax            Energy maximum, reco space
    \param r_ctbins          Number of cos-theta bins, reco space
    \param r_ctmin           cos-theta minimum, reco space
    \param r_ctmax           cos-theta maximum, reco space
    \param r_bybins          Number of bjorken-y bins, reco space
    \param r_bymin           bjorken-y minimum, reco space
    \param r_bymax           bjorken-y maximum, reco space

 */
AbsResponse::AbsResponse(reco reco_type, TString resp_name,
			 Int_t t_ebins , Double_t t_emin , Double_t t_emax ,
			 Int_t t_ctbins, Double_t t_ctmin, Double_t t_ctmax,
			 Int_t t_bybins, Double_t t_bymin, Double_t t_bymax,
			 Int_t r_ebins , Double_t r_emin , Double_t r_emax ,
			 Int_t r_ctbins, Double_t r_ctmin, Double_t r_ctmax,
			 Int_t r_bybins, Double_t r_bymin, Double_t r_bymax ) : EventFilter(reco_type) {

  fRespName = resp_name;

  vector<Double_t> t_e_edges  = NMHUtils::GetLogBins(t_ebins , t_emin , t_emax);  
  vector<Double_t> t_ct_edges = NMHUtils::GetBins   (t_ctbins, t_ctmin, t_ctmax);
  vector<Double_t> t_by_edges = NMHUtils::GetBins   (t_bybins, t_bymin, t_bymax);

  TH3D* h_bins_true = new TH3D("hbt","hbt", t_ebins , &t_e_edges[0], t_ctbins, &t_ct_edges[0], t_bybins, &t_by_edges[0]);
  
  vector<Double_t> r_e_edges  = NMHUtils::GetLogBins(r_ebins , r_emin , r_emax);  
  vector<Double_t> r_ct_edges = NMHUtils::GetBins   (r_ctbins, r_ctmin, r_ctmax);
  vector<Double_t> r_by_edges = NMHUtils::GetBins   (r_bybins, r_bymin, r_bymax);

  TH3D* h_bins_reco = new TH3D("hbr","hbr", r_ebins , &r_e_edges[0], r_ctbins, &r_ct_edges[0], r_bybins, &r_by_edges[0]);

  InitHists(fRespName, h_bins_true, h_bins_reco);

  delete h_bins_true;
  delete h_bins_reco;
}

//=========================================================================================================

/** Constructor with different binning settings for true and reco space.
    
    \param reco_type         Type of reco variables used to fill the response, e.g. DetResponse::track; See `EventFilter`
    \param resp_name         Name of the response
    \param h_bins_true       Pointer to a TH3 with desired true bins, energy on the x-axis, cos-theta on the y-axis, bjorken-y on the z-axis
    \param h_bins_reco       Pointer to a TH3 with desired reco bins, energy on the x-axis, cos-theta on the y-axis, bjorken-y on the z-axis
 */
AbsResponse::AbsResponse(reco reco_type, TString resp_name, TH3* h_bins_true, TH3* h_bins_reco) : EventFilter(reco_type) {

  fRespName = resp_name;
  InitHists(fRespName, h_bins_true, h_bins_reco);
  
}

//=========================================================================================================

/** Copy constructor 
    \param name  Name for the copy. This should not match the name of other, otherwise ROOT might get into difficulty with naming histograms.
    \param other Other instance of `AbsResponse`
*/
AbsResponse::AbsResponse(TString name, const AbsResponse &other) : EventFilter(other) {

  fRespName = name;
  InitHists(fRespName, other.fhBinsTrue, other.fhBinsReco);
  
}

//=========================================================================================================

/** Destructor */
AbsResponse::~AbsResponse() {

  if (fhBinsTrue) delete fhBinsTrue;
  if (fhBinsReco) delete fhBinsReco;
  
}

//=========================================================================================================

/** Private method to init member histograms
    \param resp_name    Name of the response, appended to histogram names
    \param h_bins_true  histogram with binning settings for true space
    \param h_bins_reco  histogram with binning settings for reco space
*/
void AbsResponse::InitHists(TString resp_name, TH3* h_bins_true, TH3* h_bins_reco) {

  fhBinsTrue = (TH3D*)h_bins_true->Clone("hbinstrue_" + resp_name);
  fhBinsTrue->SetNameTitle("hbinstrue_" + resp_name, "hbinstrue_" + resp_name);
  fhBinsTrue->SetDirectory(0);
  fhBinsTrue->Reset();

  fhBinsReco = (TH3D*)h_bins_reco->Clone("hbinsreco_" + resp_name);
  fhBinsReco->SetNameTitle("hbinsreco_" + resp_name, "hbinsreco_" + resp_name);
  fhBinsReco->SetDirectory(0);
  fhBinsReco->Reset();

}
