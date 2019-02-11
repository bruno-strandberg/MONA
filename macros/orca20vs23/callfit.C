#include "AsimovFit.h"

#include "Jeep/JParser.hh"
#include "Jeep/JMessage.hh"

int main(const int argc, const char** argv) {

  TString orca20 = "ORCA20";
  TString orca23 = "ORCA23";

  TString detector;
  Double_t th23 = 0.38;
  Bool_t InvertedOrderingData;
  TString outname;
  
  try {

    JParser<> zap("This function creates ORCA expectation value data in specified ordering and fits it assuming the other ordering; the fit results are stored in the output file. The application makes use of the AsimovFit class.");

    zap['d'] = make_field(detector, "Detector configuration, either ORCA20 or ORCA23") = std::vector< TString > { orca20, orca23 };
    zap['t'] = make_field(th23, "SinsqTh23 value at which the data is created; other osc parameters are at the central values of the specified ordering") = 0.5;
    zap['i'] = make_field(InvertedOrderingData, "Data is created with inverted ordering and fitted assuming normal ordering");
    zap['o'] = make_field(outname , "Output file where fit data is written") = "callfit_out.root";

    if ( zap.read(argc, argv) != 0 ) return 1;
  }
  catch(const exception &error) {
    FATAL(error.what() << endl);
  }

  AsimovFit::Detector det;
  if      (detector == orca23) det = AsimovFit::ORCA23;
  else if (detector == orca20) det = AsimovFit::ORCA20;
  else {
    throw std::invalid_argument("ERROR! Unknown detector type.");
  }
  
  AsimovFit AF(det);

  fitpacket data = AF.CreateData(InvertedOrderingData, th23);
  AF.FitData(data);
  AF.WriteToFile(outname, data);
  
}
