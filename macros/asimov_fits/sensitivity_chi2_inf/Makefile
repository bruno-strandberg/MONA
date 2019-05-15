
include $(MONADIR)/Makerules

CFLAGS = $(CXXFLAGS) $(ROOTCFLAGS) $(JPPCFLAGS) $(AACFLAGS)
LFLAGS = $(ROOTLIBS) $(ROOFITLIBS) $(NMHLIB) $(FITLIB)

# variable definition
#====================

# shared sources that are required by all applications
SSOURCE   = AsimovFitNBinsNOTh23Range_PercentageOfMC.C AsimovFitNBinsIOTh23Range_PercentageOfMC.C ../HelperFunctions.h

# sources of applications
PSOURCES  = AsimovFitNBinsNOTh23Range_PercentageOfMC_QSub.cpp AsimovFitNBinsIOTh23Range_PercentageOfMC_QSub.cpp
PROGRAMS  = $(PSOURCES:.cpp=)

# compilation
#=============

all: $(PROGRAMS)

$(PROGRAMS): % : %.cpp $(SSOURCE)
	@echo "  Building $@ with pre-requisities $^..."
	@g++ -O3 -o $@ $^ $(LFLAGS) $(CFLAGS)

# clean-up and print command
#===========================

.PHONY: clean all print

clean:
	@echo "  Cleaning..."
	@rm -f $(PROGRAMS)

print:
	@echo $(SSOURCE)
	@echo $(PSOURCES)
	@echo $(CFLAGS)
	@echo $(LFLAGS)
