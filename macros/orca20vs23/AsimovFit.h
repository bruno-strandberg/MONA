#ifndef AsimovFit_h
#define AsimovFit_h

#include "DetResponse.h"
#include "FitUtil.h"
#include "FitPDF.h"

#include "RooSimultaneous.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooConstVar.h"
#include "RooChi2Var.h"
#include "RooAbsReal.h"

#include "TH1.h"
#include "TCanvas.h"

//**********************************************************************************************************

/** structure to feed external chi2 terms to `MyChi2Var`; <0> is the pointer to the variable, <1> is the mean of the external measurement and <2> is the sigma of the external measurement */
typedef std::tuple<RooRealVar*, Double_t, Double_t> X2_t;

//**********************************************************************************************************

/** Dirty hack class to be able to add external chi^2 terms when performing chi^2 minimization with RooFit.

    RooFit works well with external constraints when using likelihood minimization, but the external terms are ignored when a chi2 fit is performed. This class helps to circumvent this problem.
*/
class MyChi2Var : public RooChi2Var {

 private:

  RooChi2Var *fChi2;     //!< pointer to the input `RooChi2Var`
  vector<X2_t> fExtChi2; //!< vector with data structures for external constraints

  /** This funcion returns the summed chi^2 of all of the external constraint terms in member `fExtChi2`
      \return summed chi2 of the external contraint terms
   */
  Double_t GetExternalX2() const {

    Double_t extchi2 = 0;

    for (auto x: fExtChi2) {
      Double_t data  = std::get<0>(x)->getVal();
      Double_t mean  = std::get<1>(x);
      Double_t sigma = std::get<2>(x);
      extchi2 += TMath::Power( (data - mean)/sigma, 2 );
    }

    return extchi2;
  }

 public:

  /** Default constructor
      \param chi2      A `RooChi2Var` between the model and the data
      \paran extchi2   Vector with external chi^2 term structures
   */
  MyChi2Var(RooChi2Var &chi2, vector<X2_t> &extchi2): RooChi2Var(chi2) {
    fChi2 = &chi2;       // save a pointer to the input `RooChi2Var`, need this in `evaluate()`
    fExtChi2 = extchi2;
  };

  /** Copy constructor */
  MyChi2Var(const MyChi2Var &other) : RooChi2Var( (RooChi2Var)other ) {
    fChi2 = other.fChi2;
    fExtChi2 = other.fExtChi2;
  }

  ~MyChi2Var() {};

  /** Overloaded function that adds the external chi2 terms to the input `RooChi2Var`
      Note the use of member `fChi2->getVal()`, this returns the chi^2 of the input `RooChi2Var` without external constraints. Then the external constraint terms are added manually and the sum is returned. The return of this function is the number visible to the minimizer class, e.g. `RooMinimizer`.
   */
  virtual Double_t evaluate ( ) const {
    Double_t extx2 = GetExternalX2();
    Double_t datax2 = fChi2->getVal();
    //cout << "MyChi2Var::evaluate() dataX2 and externalX2: " << datax2 << "\t" << extx2 << endl;
    return ( datax2 + extx2 );
  }

};

//**********************************************************************************************************

/** This structure helps to pass fitting data, starting values and results through various functions of the AsimovFit class and provides IO*/
struct fitpacket: public TObject {

  TString       fDetString;      //!< detector string to identify which fit this belongs to
  TH3D         *fTrkH;           //!< track expectation value
  TH3D         *fShwH;           //!< shower expectation value
  RooDataHist  *fCombH;          //!< combined dataset in RooFit
  RooArgSet    *fParData;        //!< parameter values at which the expectation data was created
  RooFitResult *fRes_1q;         //!< fit result assuming other ordering, theta-23 first quadrant
  RooFitResult *fRes_2q;         //!< fit result assuming other ordering, theta-23 second quadrant

  fitpacket(): fDetString(""), fTrkH(0), fShwH(0), fCombH(0), fParData(0), fRes_1q(0), fRes_2q(0) {};
  ~fitpacket() {};

  /** Copy constructor to help with storing data read from input on heap*/
 fitpacket(const fitpacket &other): TObject((TObject)other) {
    fDetString = other.fDetString;
    fTrkH      = (TH3D*)other.fTrkH->Clone();
    fShwH      = (TH3D*)other.fShwH->Clone();
    fCombH     = (RooDataHist*)other.fCombH->Clone();
    fParData   = (RooArgSet*)other.fParData->snapshot();
    fRes_1q    = (RooFitResult*)other.fRes_1q->Clone();
    fRes_2q    = (RooFitResult*)other.fRes_2q->Clone();

    vector< TH1* > mems = {fTrkH, fShwH};
    for (auto m: mems) m->SetDirectory(0);
  }

  /**
     Stream operator for cout.
  */
  friend std::ostream &operator << ( std::ostream &output, const fitpacket &fb ) { 
    output << "Detector: " << fb.fDetString << " Pointers: " << fb.fTrkH << " " << fb.fShwH << " " << fb.fCombH << " " << fb.fParData << " " << fb.fRes_1q << " " << fb.fRes_1q;
    return output;
  }
    
  ClassDef(fitpacket, 1)
};

/** This class and header define the functions and variables to create AsimovFit sensitvity curves in theta-23 for ORCA 23 and ORCA 20 detector geometries. */

class AsimovFit {

 public:
  
  //---------------------------------------------------------------------------
  // two supported detector types with hardcoded associated simulation and effective mass files;
  // mapping of input simulation files and effective mass files to detector type
  //---------------------------------------------------------------------------

  enum Detector { ORCA20, ORCA23 };

  std::map<Detector, TString> fDetStrings = { {ORCA20, "ORCA20"}, {ORCA23, "ORCA23"} };

  std::map<Detector, TString> fSimFiles = { 
    { ORCA20, (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_ORCA115_20x9m_ECAP1218.root" }, 
    { ORCA23, (TString)getenv("NMHDIR") + "/data/ORCA_MC_summary_ORCA115_23x9m_ECAP0418.root" } };

  std::map<Detector, TString> fEffmFiles = { 
    { ORCA20, (TString)getenv("NMHDIR") + "/data/eff_mass/EffMass_ORCA115_20x9m_ECAP1218.root" }, 
    { ORCA23, (TString)getenv("NMHDIR") + "/data/eff_mass/EffMass_ORCA115_23x9m_ECAP0418.root" } };

  //---------------------------------------------------------------------------
  // member variables
  //---------------------------------------------------------------------------
  
  // members that are populated depending on input
  TString fDetStr;
  TString fDataFile;
  TString fEffmFile;

  // responses for track and shower and relevant binning
  DetResponse *fTrkResp;
  DetResponse *fShwResp;
  
  Double_t fNebins  = 24;    // number of energy bins in range 1-100
  Double_t fNctbins = 48;    // number of cos-theta bins in range -1-1
  Double_t fNbybins = 1;     // number of bjorken-y bins in range 0-1

  // fit range and run time
  Double_t fOpTime   =  3;
  Double_t fFitEmin  =  3;
  Double_t fFitEmax  = 80;
  Double_t fFitCTmin = -1;
  Double_t fFitCTmax =  0;
  Double_t fFitBYmin =  0;
  Double_t fFitBYmax =  1;

  // data categories and members for fitting
  TString fTrkCateg = "tracks";
  TString fShwCateg = "showers";

  RooCategory     *fCategs;
  FitUtil         *fFitUtil;
  FitPDF          *fTrkPdf;
  FitPDF          *fShwPdf;
  RooSimultaneous *fSimPdf; // a joint pdf for tracks and showers, does not include external constraints
  RooProdPdf      *fFitPdf; // pdf that includes external constraints, useful for `LikelihoodScan` and `Contour`

  // variables for an external constraint on th13
  RooGaussian     *fTh13C;
  static constexpr double fTh13mean  = 0.0215;   // from PDG
  static constexpr double fTh13sigma = 0.0025/3; // 3sigma limits are 0.019, 0.024, i.e. += 0.0025
  vector<X2_t> fExternalX2terms;                 // vector with external chi^2 terms
  
  // vector for storing fitpacket's read from a file
  vector< fitpacket* > fFPs;
  
  static constexpr double fNOTSET = 1000.; // variable to indicate un-used function parameters

  //---------------------------------------------------------------------------
  // functions
  //---------------------------------------------------------------------------

  AsimovFit(Detector detector);
  ~AsimovFit();

  void ParseInput(Detector detector);
  void InitResponses(Detector detector);
  void InitPdfs();

  fitpacket CreateData(Bool_t NOdata = kTRUE, 
				  Double_t th23 = fNOTSET, Double_t th12 = fNOTSET, Double_t th13 = fNOTSET, 
				  Double_t dcp = fNOTSET, Double_t dm21 = fNOTSET, Double_t dm31 = fNOTSET);
  void FitData(fitpacket& fp);
  void WriteToFile(TString outfile, vector<fitpacket> fps);
  void WriteToFile(TString outfile, fitpacket fp) { WriteToFile(outfile, vector<fitpacket> {fp} ); }
  void ReadFromFile(TString filename);
  void LikelihoodScan(TString var_name = "SinsqTh23", Double_t var_min = 0.35, 
		      Double_t var_max = 0.65, Double_t th23 = 0.6);
  void Contour(Double_t th23_data, TString var1="SinsqTh23", TString var2 = "SinsqTh13");
  RooDataHist* CombineData(TH1* trkH, TH1* shwH);
  void SetModelToFitresult(RooFitResult *fr);
  std::tuple<Double_t, Double_t, Double_t> GetChi2(fitpacket &fp, Bool_t q1);
  std::tuple<Double_t, Double_t, Double_t> GetAsym(fitpacket &fp, Bool_t q1);

  fitpacket* FindPacket(Double_t sinsqth23, Detector det);

};

#endif
