#ifndef HELPERS_ASIMOV
#define HELPERS_ASIMOV

#include "TH3.h"
#include "TMath.h"

#include "FitUtil.h"


/**
 *  Function calculate the chi2-value between two TH1. 
 *  
 *
 * \param h1           Histogram
 * \param h2           Histogram
 * \return             Chi2 value between two histograms.
 */
Double_t HistoChi2Test(TH1* h1, TH1* h2,
                       Double_t xlow = -1e10, Double_t xhigh = 1e10,
                       Double_t ylow = -1e10, Double_t yhigh = 1e10) {

  Double_t chi  = 0;
  Double_t chi2 = 0;
  if (h1->GetDimension() != h2->GetDimension()) {
    throw std::invalid_argument( "ERROR! input histograms dimensions mismatch.");
  }

  for (Int_t xb = 1; xb <= h1->GetXaxis()->GetNbins(); xb++) {
    for (Int_t yb = 1; yb <= h1->GetYaxis()->GetNbins(); yb++) {
      for (Int_t zb = 1; zb <= h1->GetZaxis()->GetNbins(); zb++) {
        Double_t val1 = h1->GetBinContent(xb, yb, zb);
        Double_t val2 = h2->GetBinContent(xb, yb, zb);
        
        if (val1 == 0) {
          chi = 0;
        } else {
          chi = (val1 - val2) / TMath::Sqrt(val1) ;
        }
        Double_t xc = h1->GetXaxis()->GetBinCenter(xb);
        Double_t yc = h1->GetYaxis()->GetBinCenter(yb);

        if ( (xc < xlow) || (xc > xhigh) || ( yc < ylow) || ( yc > yhigh ) ) {
          chi = 0.;
        }
        chi2 += chi*chi;
      }
    }
  }

  return TMath::Sqrt(chi2);
}


Double_t HistoChi2Calc(TH1* h1, TH1* h2,
                       Double_t xlow = -1e10, Double_t xhigh = 1e10,
                       Double_t ylow = -1e10, Double_t yhigh = 1e10,
                       Double_t zlow = -1e10, Double_t zhigh = 1e10) {

  if (h1->GetDimension() != h2->GetDimension()) {
    throw std::invalid_argument( "ERROR! input histograms dimensions mismatch.");
  }

  TH1* _h1 = (TH1*)h1->Clone("h1clone");
  TH1* _h2 = (TH1*)h2->Clone("h2clone");
  _h1->SetDirectory(0);
  _h2->SetDirectory(0);

  for (Int_t xb = 1; xb <= h1->GetXaxis()->GetNbins(); xb++) {
    for (Int_t yb = 1; yb <= h1->GetYaxis()->GetNbins(); yb++) {
      for (Int_t zb = 1; zb <= h1->GetZaxis()->GetNbins(); zb++) {
        Double_t xbincenter = _h1->GetXaxis()->GetBinCenter(xb);
        Double_t ybincenter = _h1->GetYaxis()->GetBinCenter(yb);
        if ((xbincenter < xlow) or (ybincenter < ylow)  or
            (xbincenter > xhigh) or (ybincenter > yhigh)) {

          _h1->SetBinContent(xb, yb, zb, 0);
          _h1->SetBinError(xb, yb, zb, 0);
          _h2->SetBinContent(xb, yb, zb, 0);
          _h2->SetBinError(xb, yb, zb, 0);
        }
      }
    }
  }

  return _h1->Chi2Test(_h2, "WW CHI2");
}


// set up a pointer to a function to call resetting the parameters to central values
void ResetToCentral(FitUtil& fitutil) {
  RooRealVar *dm31 = (RooRealVar*)fitutil.GetVar("Dm31");
  if ( !dm31 ) throw std::logic_error("ERROR! ResetToCentral cannot find dm31");

  if ( dm31->getVal() > 0 ) { 
    cout << "NOTICE: Set values to IO " << endl; 
    fitutil.SetIOcentvals(); 
  }
  else { 
    cout << "NOTICE: Set values to NO " << endl; 
    fitutil.SetNOcentvals(); 
  }
}



/** Function to set the oscillation parameter limits corresponding to normal mass ordering with 
 *  \f$ sin^2\theta_{12} \f$ value to free s.t. you can fit in both quantiles.
 */
void SetNOlimsChi2Fit(FitUtil* fitutil) {

  fitutil->SetNOlims();

  // If you want to fit both quantiles of th23 you need free limits here
  fitutil->GetVar("SinsqTh23")->setMin(0.);
  fitutil->GetVar("SinsqTh23")->setMax(1.);
}


/** Function to set the oscillation parameter limits corresponding to inverted mass ordering with 
 *  \f$ sin^2\theta_{12} \f$ value to free s.t. you can fit in both quantiles.
 */
void SetIOlimsChi2Fit(FitUtil* fitutil) {

  fitutil->SetIOlims();

  // If you want to fit both quantiles of th23 you need free limits here
  fitutil->GetVar("SinsqTh23")->setMin(0.);
  fitutil->GetVar("SinsqTh23")->setMax(1.);
}


/** Function to get the folders to save the detector responses and root files
 */
TString DetectorResponseFolder(Int_t n) {
  TString filefolder = Form("./detector_responses/pid_binning_%i/", n);
  cout << "NOTICE: Saving into folder " << filefolder << endl;
  return filefolder;
}


/** Function to create a map of the PID bin edges
 */
std::map<Int_t, Double_t> SetPIDCase(Int_t n) {
  std::map<Int_t, Double_t> pid_map;

  switch (n) {
    case 2:
      cout << "NOTICE: Set PID case " << n << endl;
      pid_map.insert(std::make_pair(0, 0.0)); // shower
      pid_map.insert(std::make_pair(1, 0.6)); // track
      pid_map.insert(std::make_pair(2, 1.0)); // upper limit
      break;
    case 3:
      cout << "NOTICE: Set PID case " << n << endl;
      pid_map.insert(std::make_pair(0, 0.0)); // shower
      pid_map.insert(std::make_pair(1, 0.4)); // middle group: shower
      pid_map.insert(std::make_pair(2, 0.6)); // track
      pid_map.insert(std::make_pair(3, 1.0)); // upper limit 
      break;
    case 4:
      cout << "NOTICE: Set PID case " << n << endl;
      pid_map.insert(std::make_pair(0, 0.0)); // shower
      pid_map.insert(std::make_pair(1, 0.4)); // middle group: shower
      pid_map.insert(std::make_pair(2, 0.6)); // middle group: track
      pid_map.insert(std::make_pair(3, 0.8)); // track
      pid_map.insert(std::make_pair(4, 1.0)); // upper limit 
      break;
    case 5:
      cout << "NOTICE: Set PID case " << n << endl;
      // The for loop also has to create an upper limit, thus +1
      for (Int_t i = 0; i < n+1; i++) { 
        pid_map.insert(std::make_pair(i, i / float(n) )); // This is i * PID_STEP
      }
      break;
    case 10:
      cout << "NOTICE: Set PID case " << n << endl;
      for (Int_t i = 0; i < n+1; i++) { 
        pid_map.insert(std::make_pair(i, i / float(n) )); 
      }
      break;
    case 21:
      cout << "NOTICE: Set PID case " << n << endl;
      pid_map.insert(std::make_pair(0, 0.0)); // shower
      pid_map.insert(std::make_pair(1, 0.7)); // track
      pid_map.insert(std::make_pair(2, 1.0)); // upper limit
      break;
    case 22:
      cout << "NOTICE: Set PID case " << n << endl;
      pid_map.insert(std::make_pair(0, 0.0)); // shower
      pid_map.insert(std::make_pair(1, 0.8)); // track
      pid_map.insert(std::make_pair(2, 1.0)); // upper limit
      break;
    case 31:
      cout << "NOTICE: Set PID case " << n << endl;
      cout << "NOTICE: Non-standard ranges on the PID categories!" << endl;
      cout << "NOTICE: The user must set the track-shower cut to 0.8" << endl;
      pid_map.insert(std::make_pair(0, 0.0)); // shower
      pid_map.insert(std::make_pair(1, 0.4)); // middle group: shower
      pid_map.insert(std::make_pair(2, 0.8)); // track
      pid_map.insert(std::make_pair(3, 1.0)); // upper limit 
      break;
    case 32:
      cout << "NOTICE: Set PID case " << n << endl;
      cout << "NOTICE: Non-standard ranges on the PID categories!" << endl;
      cout << "NOTICE: The user must set the track-shower cut to 0.7" << endl;
      pid_map.insert(std::make_pair(0, 0.0)); // shower
      pid_map.insert(std::make_pair(1, 0.3)); // middle group: shower
      pid_map.insert(std::make_pair(2, 0.7)); // track
      pid_map.insert(std::make_pair(3, 1.0)); // upper limit 
      break;
    case 41:
      cout << "NOTICE: Set PID case " << n << endl;
      cout << "NOTICE: Non-standard ranges on the PID categories!" << endl;
      cout << "NOTICE: The user must set the track-shower cut to 0.7" << endl;
      pid_map.insert(std::make_pair(0, 0.0)); // shower
      pid_map.insert(std::make_pair(1, 0.3)); // lower group: shower
      pid_map.insert(std::make_pair(2, 0.5)); // middle group: shower
      pid_map.insert(std::make_pair(3, 0.7)); // track
      pid_map.insert(std::make_pair(4, 1.0)); // upper limit 
      break;
  }

  return pid_map;

}


/** Function to give the PID cut for track and shower.
 *  
 *  There are some exceptional cases for the PID binnins, they are stored under
 *  the values 31, 42, etc. These number mean: 
 *  31 --> 3 PID bins, case 1
 *  42 --> 4 PID bins, case 2, etc.
 *  This naming scheme works under the assumption that we never take more than 10 
 *  PID categories.
 *
 *  returns Int_t below 10.
 */
Int_t GetNumPIDCats(Int_t n) {
  if (n <= 10) return n;
  else return n / 10;
}

/** Funtion to return the PID cut between tracks and showers.
 *  Watch out that these values are hardcoded and depend on the
 *  values above. Above is safer since it does not depend on
 *  anything. This function depends on SetPIDCase.
 */
Double_t GetPIDCut(Int_t n) {
  if (n <= 10) return 0.6;
  else if (n == 21) return 0.7;
  else if (n == 22) return 0.8;
  else if (n == 31) return 0.8;
  else if (n == 32) return 0.7;
  else if (n == 41) return 0.7;
  else {
    cout << "ERROR: Invalid PID number given. Exiting." << endl;   
    exit(-1);
  }
}


/** Function to create a map of the PID bin edges
 */
std::vector< std::tuple<Double_t, Double_t> > GetEnergyRanges(Int_t n) {
  std::vector< std::tuple<Double_t, Double_t> > range_map;

  switch (n) {
    case 2:
      cout << "NOTICE: Energy range case " << n << endl;
      // The values of all these maps are based on 10% MC Uncertainty as cut-off for NOT using the
      // events by using visual inspection. So outside these ranges, the MCU on the events is >10%
      //
      // =================================================
      //  THESE VALUES HAVE TO BE RECHECKED AS OF 20190403
      // =================================================
      //
      range_map.push_back(std::make_tuple(2.0, 80.0)); // shower
      range_map.push_back(std::make_tuple(3.0, 80.0)); // track
      break;
    case 3:
      cout << "NOTICE: Energy range case " << n << endl;
      range_map.push_back(std::make_tuple(2.0, 80.0)); // shower
      range_map.push_back(std::make_tuple(2.0, 30.0)); // middle group: shower
      range_map.push_back(std::make_tuple(3.0, 80.0)); // track
      break;
    case 4:
      cout << "NOTICE: Energy range case " << n << endl;
      range_map.push_back(std::make_tuple(2.0, 80.0)); // shower
      range_map.push_back(std::make_tuple(2.0, 30.0)); // middle group: shower
      range_map.push_back(std::make_tuple(3.0, 15.0)); // middle group: track (kind of)
      range_map.push_back(std::make_tuple(5.0, 80.0)); // track
      break;
    case 5:
      cout << "NOTICE: Energy range case " << n << endl;
      range_map.push_back(std::make_tuple(2.0, 80.0)); // shower
      range_map.push_back(std::make_tuple(2.0, 50.0));
      range_map.push_back(std::make_tuple(2.0, 30.0));
      range_map.push_back(std::make_tuple(3.0, 15.0)); // track
      range_map.push_back(std::make_tuple(5.0, 80.0));
      break;
    case 10:
      cout << "NOTICE: Energy range case " << n << endl;
      range_map.push_back(std::make_tuple(7.0, 80.0)); // shower
      range_map.push_back(std::make_tuple(2.5, 70.0));
      range_map.push_back(std::make_tuple(2.5, 50.0));
      range_map.push_back(std::make_tuple(2.0, 30.0));
      range_map.push_back(std::make_tuple(2.0, 20.0));
      range_map.push_back(std::make_tuple(2.0, 15.0));
      range_map.push_back(std::make_tuple(3.0, 15.0)); // track
      range_map.push_back(std::make_tuple(3.5, 15.0));
      range_map.push_back(std::make_tuple(5.0, 15.0));
      range_map.push_back(std::make_tuple(5.0, 80.0));
      break;
  }

  return range_map;

}


/**
 *  Function calculate the `difference` between two TH1. 
 *  
 *
 * \param h1           Histogram
 * \param h2           Histogram
 * \return             Histogram containing difference between two histograms. Difference is 
 *                     defined as h1 - h2 / sigma_h1
 */
TH3D* HistoDifference(TH1* h1, TH1* h2) {

  TH3D* ret = (TH3D*)h1->Clone("diff_histo");
  ret->Reset();
  ret->SetDirectory(0);

  if (h1->GetDimension() != h2->GetDimension()) {
    throw std::invalid_argument( "ERROR! input histograms dimensions mismatch.");
  }

  for (Int_t xb = 1; xb <= h1->GetXaxis()->GetNbins(); xb++) {
    for (Int_t yb = 1; yb <= h1->GetYaxis()->GetNbins(); yb++) {
      for (Int_t zb = 1; zb <= h1->GetZaxis()->GetNbins(); zb++) {
        Double_t binval;
        Double_t val1 = h1->GetBinContent(xb, yb, zb);
        Double_t val2 = h2->GetBinContent(xb, yb, zb);
        Double_t sig1 = h1->GetBinError(xb, yb, zb);
        
        if (sig1 == 0) {
          binval = 0;
        } else {
          binval = (val1 - val2) / sig1 ;
        }

        ret->SetBinContent(xb, yb, zb, binval);

      }
    }
  }

  return ret;
}

#endif
