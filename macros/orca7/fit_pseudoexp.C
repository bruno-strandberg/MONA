#include "NMHUtils.h"
#include "ORCA7.h"
#include "ORCA7.C"

using namespace RooFit;

//***************************************************************************************

void FixSystematics(FitUtil *fu) {

  fu->GetVar("E_tilt")     ->setConstant( kTRUE );
  fu->GetVar("ct_tilt")    ->setConstant( kTRUE );
  fu->GetVar("skew_mu_amu")->setConstant( kTRUE );
  fu->GetVar("skew_e_ae")  ->setConstant( kTRUE );
  fu->GetVar("skew_mu_e")  ->setConstant( kTRUE );
  fu->GetVar("NC_norm")    ->setConstant( kTRUE );
  fu->GetVar("Tau_norm")   ->setConstant( kTRUE );
  fu->GetVar("E_scale")    ->setConstant( kTRUE );
 
}

//***************************************************************************************

void SystematicsToCentVals(FitUtil *fu) {

  fu->GetVar("E_tilt")     ->setVal( 0. );
  fu->GetVar("ct_tilt")    ->setVal( 0. );
  fu->GetVar("skew_mu_amu")->setVal( 0. );
  fu->GetVar("skew_e_ae")  ->setVal( 0. );
  fu->GetVar("skew_mu_e")  ->setVal( 0. );
  fu->GetVar("NC_norm")    ->setVal( 1. );
  fu->GetVar("Tau_norm")   ->setVal( 1. );
  fu->GetVar("E_scale")    ->setVal( 0. );

}

//***************************************************************************************

/* 
   Info here
*/
void fit_pseudoexp(TString output        = "fit_pseudoexp.root",
		   Double_t sinsqth23    = 0.425,
		   Bool_t fixSystematics = kTRUE,
		   Int_t ncpu            = 1,
		   Int_t seed            = 418) {

  //===============================================================================
  // set up pdf's
  //===============================================================================

  ORCA7 o7(kTRUE);

  FitUtilWsyst futil = *o7.fFitUtil;
  std::map<TString, FitPDF*> pdfs = o7.fPdfs;
  for (auto &P: pdfs) { P.second->SetSeed(seed); } // seed of pseudo-experiment generation

  //===============================================================================
  // create pseudo-experiments
  //===============================================================================

  futil.SetNOcentvals();
  futil.GetVar("SinsqTh23")->setVal( sinsqth23 );

  std::map<string, TH1*> exph;
    
  for (auto P: pdfs) {
    string  key = (string)P.first;
    FitPDF *pdf = P.second;
    TH3D   *exp = pdf->SimplePseudoExp(TString("h_") + P.first, kFALSE);
    exph.insert( std::make_pair(key, (TH1*)exp) );
  }

  RooArgSet *datapars = (RooArgSet*)futil.GetSet().snapshot(kTRUE);

  //===============================================================================
  // configure simultaneous fitting
  //===============================================================================

  RooCategory categ("categ","categ");
  for (auto exp: exph) { categ.defineType( (TString)exp.first ); }

  RooDataHist combh("combdata","combined data", futil.GetObs(), categ, exph);

  RooSimultaneous simPdf("simPdf","simultaneous pdf", categ);
  for (auto &P: pdfs) { simPdf.addPdf( *P.second, P.first ); }

  // need to add some priors; skew parameter priors ball-parked from Barr
  // the uncertainty there is <=20%, a prior of 20% hopefully helps to avoid the fitter going
  // to regions we know are not physical
  RooGaussian mu_amu_prior("mu_amu_prior", "mu_amu_prior", *futil.GetVar("skew_mu_amu"), RooConst(0.), RooConst(0.2) );
  RooGaussian e_ae_prior  ("e_ae_prior"  , "e_ae_prior"  , *futil.GetVar("skew_e_ae")  , RooConst(0.), RooConst(0.2) );
  RooGaussian mu_e_prior  ("mu_e_prior"  , "mu_e_prior"  , *futil.GetVar("skew_mu_e")  , RooConst(0.), RooConst(0.2) );

  // energy scale prior is a complete guess; setting 0.15, which means 2sigma is 30%, seems kind-of reasonable
  RooGaussian escale_prior("escale_prior", "escale_prior", *futil.GetVar("E_scale"), RooConst(0.), RooConst(0.15));

  // nc normalisation taken from Neutrino2018 poster
  RooGaussian ncnorm_prior("ncnorm_prior", "ncnorm_prior", *futil.GetVar("NC_norm"), RooConst(0.), RooConst(0.1) );

  RooArgSet components( simPdf, mu_amu_prior, e_ae_prior, mu_e_prior, escale_prior, ncnorm_prior );
  RooProdPdf fitPdf("fitPdf","fitPdf", components);

  //===============================================================================
  // fix some parameters, create th23 ranges
  //===============================================================================

  futil.GetVar("SinsqTh12")->setConstant(kTRUE);
  futil.GetVar("SinsqTh13")->setConstant(kTRUE);
  futil.GetVar("dcp")->setConstant(kTRUE);
  futil.GetVar("Dm21")->setConstant(kTRUE);
  if (fixSystematics) FixSystematics(&futil);

  // fix some additional parameters
  futil.GetVar("Tau_norm")->setConstant(kTRUE);
  futil.GetVar("ct_tilt")->setConstant(kTRUE);

  futil.GetVar("SinsqTh23")->setRange("firstq" , 0. , 0.5  );
  futil.GetVar("SinsqTh23")->setRange("secondq", 0.5, 1.   );

  //===============================================================================
  // fit in first and second quadrant
  //===============================================================================

  // fitting in first quadrant
  futil.SetNOcentvals();
  SystematicsToCentVals(&futil);
  futil.GetVar("SinsqTh23")->setVal(0.4);
  RooFitResult *fitres_q1 = fitPdf.fitTo( combh, Save(kTRUE), SumW2Error(kFALSE), Range("firstq" ), NumCPU(ncpu) );

  // fitting in second quadrant
  futil.SetNOcentvals();
  SystematicsToCentVals(&futil);
  futil.GetVar("SinsqTh23")->setVal(0.6);
  RooFitResult *fitres_q2 = fitPdf.fitTo( combh, Save(kTRUE), SumW2Error(kFALSE), Range("secondq"), NumCPU(ncpu) );

}
