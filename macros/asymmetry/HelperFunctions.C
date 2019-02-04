#include <iostream> 
#include <fstream> 


 /* IO Functions
  */

 /* Small function to read a detector response file and take the TH3D
  * responses of tracks and showers and project them into the 2D E-ct
  * plane
  */
std::tuple<TH2D*, TH2D*, TH2D*> ReadDetectorResponseFile(TString filename) {

  TFile *f = TFile::Open(filename, "READ");

  TH3D *h3_t = (TH3D*)f->Get("detected_tracks")->Clone();
  TH3D *h3_s = (TH3D*)f->Get("detected_showers")->Clone();
  TH3D *h3_m = (TH3D*)f->Get("detected_mc")->Clone();
  h3_t->SetDirectory(0);
  h3_s->SetDirectory(0);
  h3_m->SetDirectory(0);

  TH2D *h_t = (TH2D*)h3_t->Project3D("yx");
  TH2D *h_s = (TH2D*)h3_s->Project3D("yx");
  TH2D *h_m = (TH2D*)h3_m->Project3D("yx");

  return std::make_tuple(h_t, h_s, h_m);
}


 /* Small function to read a detector response file and take the TH3D
  * responses of tracks and showers and project them into the 2D E-ct
  * plane.
  *
  * This is an old version using only track and shower. To be deprecated.
  */
std::tuple<TH2D*, TH2D*> ReadDetectorResponseFile2(TString filename) {

  TFile *f = TFile::Open(filename, "READ");

  TH3D *h3_t = (TH3D*)f->Get("detected_tracks")->Clone();
  TH3D *h3_s = (TH3D*)f->Get("detected_showers")->Clone();
  h3_t->SetDirectory(0);
  h3_s->SetDirectory(0);

  TH2D *h_t = (TH2D*)h3_t->Project3D("yx");
  TH2D *h_s = (TH2D*)h3_s->Project3D("yx");

  return std::make_tuple(h_t, h_s);
}


 /* Small function to read a correlation coefficient file and take the 
  * TH3D and project them into the 2D E-ct plane
  */
std::tuple<TH2D*, TH2D*, TH2D*> ReadCorrelationFile(TString filename) {

  TFile *f = TFile::Open(filename, "READ");

  TH3D *h3_corr_t = (TH3D*)f->Get("correlation_track_bins")->Clone();
  TH3D *h3_corr_s = (TH3D*)f->Get("correlation_shower_bins")->Clone();
  TH3D *h3_corr_m = (TH3D*)f->Get("correlation_mc_bins")->Clone();
  h3_corr_t->SetDirectory(0);
  h3_corr_s->SetDirectory(0);
  h3_corr_m->SetDirectory(0);

  TH2D *h_t = (TH2D*)h3_corr_t->Project3D("yx");
  TH2D *h_s = (TH2D*)h3_corr_s->Project3D("yx");
  TH2D *h_m = (TH2D*)h3_corr_m->Project3D("yx");

  return std::make_tuple(h_t, h_s, h_m);
}


  /* Printer functions
   */

  /* Function to print asymmetry with errors.
   */
void PrintAsymmetryWithErrors(string type, double asym, double asym_err) {
    cout << "Asymmetry for " << type << ": " << asym << " +- " << asym_err << " (" << 100*asym_err/asym << "%)" << endl;
}

void PrintAsymmetryWithErrors(string type, double asym, double asym_err, TString filename) {
    cout << "Asymmetry for " << type << ": " << asym << " +- " << asym_err << " (" << 100*asym_err/asym << "%)" << endl;

    ofstream outputfile(filename, std::ios_base::app);
    outputfile << "Asymmetry for " << type << ": " << asym << " +- " << asym_err << " (" << 100*asym_err/asym << "%)" << endl;
    outputfile.close();
}


  /* Function to print asymmetry with errors and events counts.
   */
void PrintAsymmetryWithErrors(string type, double asym, double asym_err, double count_NO, double count_IO) {
    cout << "Asymmetry for " << type << ": " << asym << " +- " << asym_err << " (" << 100*asym_err/asym << "%)" << "  [" << count_NO << " : " << count_IO << "]" << endl;
}


  /* Mathematical functions
   */

  /* Function to calculated the normed innerproduct of two vectors. 
   * This is used to calculate the correlation coefficient between
   * two vectors, which is equivalent to the innerproduct.
   */
Double_t normed_innerproduct(std::vector<Double_t> a, std::vector<Double_t> b) { // This gives the correlation coefficient between vectors
  Double_t inner  = std::inner_product(std::begin(a), std::end(a), std::begin(b), 0.0);
  Double_t norm_a = std::sqrt(std::inner_product(std::begin(a), std::end(a), std::begin(a), 0.0));
  Double_t norm_b = std::sqrt(std::inner_product(std::begin(b), std::end(b), std::begin(b), 0.0));
  if (norm_a * norm_b != 0) {
    Double_t inner_normed = inner / (norm_a * norm_b);
    return inner_normed;
  } else {
    return 0;
  }
}

  /* Function to normalize TH2D objects in the Y direction per column.
   */
void GetNormalizedSlicesY(TH2D* input) {
  for (int i_b = 0; i_b <= input->GetXaxis()->GetNbins()+1; i_b++) {
    TH1D* h_sl = input->ProjectionY("h_sl", i_b, i_b);
    h_sl->SetDirectory(0);
    double integral = h_sl->Integral();
    if (integral > 0) {
      for (int i_y = 0; i_y <= input->GetYaxis()->GetNbins()+1; i_y++) {
        if (input->GetBinContent(i_b, i_y) > 0){
          input->SetBinContent(i_b, i_y, input->GetBinContent(i_b, i_y) / integral);
        }
        else{
          input->SetBinContent(i_b, i_y, 0);
        }
      }
    }
  }
}

  /* Function to return a histogram with the errors of input histogram.
   */
TH2D* GetErrorHisto2D(TH2D* input) {
  Int_t nbinsX = input->GetXaxis()->GetNbins() + 2;
  Int_t nbinsY = input->GetYaxis()->GetNbins() + 2;
  TString name = input->GetName();
  name = name + "_err";
  TH2D* output = (TH2D*)input->Clone(name);
  for (int i_b = 0; i_b <= nbinsX * nbinsY; i_b++) {
    output->SetBinContent(i_b, input->GetBinError(i_b) );
    output->SetBinError(i_b, 0);
  }
  return output;
}

  /* Function to return a vector containing the Y slices of a TH2D
   */
std::vector<TH1D*> VectorProjectionY2D(TH2D* input) {
  std::vector<TH1D*> output;
  for(Int_t i = 1; i < input->GetYaxis()->GetNbins() + 1; i++) {
    TString proj_name = input->GetName();
    proj_name = proj_name + std::to_string(i);
    TH1D* proj = input->ProjectionY(proj_name, i, i);
    proj->SetDirectory(0);
    output.push_back(proj);
  }
  return output;
}

  /* Function to get the under- and overflow counts of a TH2D.
   */
std::tuple<Double_t, Double_t, Double_t, Double_t> GetOverFlow2D(TH2D* input) {
  TH1D* h_proj_x = input->ProjectionX("h_proj_x");
  TH1D* h_proj_y = input->ProjectionY("h_proj_y");

  Double_t x_underflow = h_proj_x->GetBinContent(0);
  Double_t x_overflow  = h_proj_x->GetBinContent(h_proj_x->GetNbinsX()+1); // It seems always to project to X.
  Double_t y_underflow = h_proj_y->GetBinContent(0);
  Double_t y_overflow  = h_proj_y->GetBinContent(h_proj_y->GetNbinsX()+1);

  return std::make_tuple(x_underflow, x_overflow, y_underflow, y_overflow);
}


  /* Function to get the under- and overflow counts of a TH2D.
   */
void GetValues3D(TH3D* input) {

  for (Int_t ebin = 1; ebin <= input->GetXaxis()->GetNbins(); ebin++) {
    for (Int_t ctbin = 1; ctbin <= input->GetYaxis()->GetNbins(); ctbin++) {
      for (Int_t bybin = 1; bybin <= input->GetZaxis()->GetNbins(); bybin++) {

        Double_t E  = input->GetXaxis()->GetBinCenter( ebin );
        Double_t ct = input->GetYaxis()->GetBinCenter( ctbin );
        Double_t by = input->GetZaxis()->GetBinCenter( bybin );

        Double_t bincontent = input->GetBinContent( ebin, ctbin, bybin );
        Double_t binerror   = input->GetBinError(   ebin, ctbin, bybin );

        if (bincontent != 0) {
          cout << "Bincontent at " << ebin << " " << ctbin << " " << bybin << " is " << bincontent << " +- " << binerror << "   " << bincontent/binerror << endl;
        }
      }
    }
  }
}

  /* Plotting functions
   */

void SetPlotTitle(TH1* h1, TString title) {
  h1->SetTitle(title);
  h1->SetXTitle("E_{reco} [GeV]");
  h1->SetYTitle("cos(#theta_{reco})");
}

void SetLabelSizes(TH1* h1, Double_t size) {
  h1->GetXaxis()->SetTitleSize(size);
  h1->GetXaxis()->SetTitleOffset(1.2);
  h1->GetXaxis()->SetTitleFont(62);
  h1->GetXaxis()->SetLabelFont(62);
  h1->GetXaxis()->SetLabelSize(size);
  h1->GetYaxis()->SetTitleSize(size);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTitleFont(62);
  h1->GetYaxis()->SetLabelFont(62);
  h1->GetYaxis()->SetLabelSize(size);
}
