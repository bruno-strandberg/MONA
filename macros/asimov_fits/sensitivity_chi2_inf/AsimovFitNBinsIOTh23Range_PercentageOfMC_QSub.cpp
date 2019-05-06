#include <iostream>
#include "AsimovQSubHeader.h"
      
int main(int argc, char* argv[]) {
 
  if (argc != 3) {
    std::cerr << "ERROR: Wrong number of arguments, must be 2: job number and PID binning" << std::endl;
    exit(-1);
  }
  Int_t jobNum = std::atoi(argv[1]);
  Int_t pidNum = std::atoi(argv[2]);
  AsimovFitNBinsIOTh23Range_PercentageOfMC(jobNum, pidNum);
  return 0;
}
