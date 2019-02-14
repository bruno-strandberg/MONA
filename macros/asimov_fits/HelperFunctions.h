#include "TH3.h"
#include "TMath.h"

#include "FitUtil.h"


//****************************************************************************

/**
 *  Function calculate the chi2-value between two TH1. 
 *  
 *
 * \param h1           Histogram
 * \param h2           Histogram
 * \return             Chi2 value between two histograms.
 */
Double_t HistoChi2Test(TH3D* h1, TH3D* h2,
                       Double_t xlow = -1e10, Double_t xhigh = 1e10,
                       Double_t ylow = -1e10, Double_t yhigh = 1e10) {

  Double_t chi  = 0;
  Double_t chi2 = 0;
  if (h1->GetDimension() != h2->GetDimension()) {
    throw std::invalid_argument( "ERROR! input histograms dimensions mismatch.");
  }

  for (Int_t xb = 1; xb <= h1->GetXaxis()->GetNbins(); xb++) {
    for (Int_t yb = 1; yb <= h1->GetYaxis()->GetNbins(); yb++) {
      for (Int_t zb = 1; zb <= h1->GetYaxis()->GetNbins(); zb++) {
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

  fitutil->GetVar("SinsqTh12")->setMin(f_NO_sinsqth12.min);
  fitutil->GetVar("SinsqTh12")->setMax(f_NO_sinsqth12.max);
  
  fitutil->GetVar("SinsqTh13")->setMin(f_NO_sinsqth13.min);
  fitutil->GetVar("SinsqTh13")->setMax(f_NO_sinsqth13.max);

  // If you want to fit both quantiles of th23 you need free limits here
  fitutil->GetVar("SinsqTh23")->setMin(0.);
  fitutil->GetVar("SinsqTh23")->setMax(1.);

  fitutil->GetVar("dcp")->setMin(f_NO_dcp.min);
  fitutil->GetVar("dcp")->setMax(f_NO_dcp.max);

  fitutil->GetVar("Dm21")->setMin(f_NO_dm21.min);
  fitutil->GetVar("Dm21")->setMax(f_NO_dm21.max);
  
  fitutil->GetVar("Dm31")->setMax(f_NO_dm31.max);
  fitutil->GetVar("Dm31")->setMin(f_NO_dm31.min);
}


/** Function to set the oscillation parameter limits corresponding to inverted mass ordering with 
 *  \f$ sin^2\theta_{12} \f$ value to free s.t. you can fit in both quantiles.
 */
void SetIOlimsChi2Fit(FitUtil* fitutil) {

  fitutil->GetVar("SinsqTh12")->setMin(f_IO_sinsqth12.min);
  fitutil->GetVar("SinsqTh12")->setMax(f_IO_sinsqth12.max);
  
  fitutil->GetVar("SinsqTh13")->setMin(f_IO_sinsqth13.min);
  fitutil->GetVar("SinsqTh13")->setMax(f_IO_sinsqth13.max);

  // If you want to fit both quantiles of th23 you need free limits here
  fitutil->GetVar("SinsqTh23")->setMin(0.);
  fitutil->GetVar("SinsqTh23")->setMax(1.);

  fitutil->GetVar("dcp")->setMin(f_IO_dcp.min);
  fitutil->GetVar("dcp")->setMax(f_IO_dcp.max);

  fitutil->GetVar("Dm21")->setMin(f_IO_dm21.min);
  fitutil->GetVar("Dm21")->setMax(f_IO_dm21.max);
  
  fitutil->GetVar("Dm31")->setMin(f_IO_dm31.min);
  fitutil->GetVar("Dm31")->setMax(f_IO_dm31.max);
}
