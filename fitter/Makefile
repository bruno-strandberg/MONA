
#===============================================================
# setup cxx and library flags
#===============================================================
JPPCFLAGS     = -I$(JPP_DIR)/software
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lRooFit -lRooFitCore
NMHLIBS       = -L $(NMHDIR)/common_software -lnmhsoft -L $(OSCPROBDIR) -lOscProb

CXX           = g++
CXXFLAGS      = -g -Wall -std=c++11 -fPIC $(ROOTCFLAGS) $(JPPCFLAGS)
LDFLAGS       = $(ROOTLIBS) $(NMHLIBS)

#===============================================================
# define variables for complilation
#===============================================================
PACKAGES = Fitter Asymmetry
COMMON_HEADERS = $(filter-out $(PACKAGES:=.h), $(wildcard *.h))
COMMON_SOURCES = $(filter-out $(PACKAGES:=.C), $(wildcard *.C))
DICTIONARY = dict.cxx

#===============================================================
# compilation
#===============================================================

all: $(PACKAGES)

$(PACKAGES): %: %.C $(DICTIONARY) $(COMMON_SOURCES)
	@echo "  Building $@ with pre-requisities $^..."
	@g++ -O3 -o $@ $^ $(LDFLAGS) $(CXXFLAGS)

$(DICTIONARY): $(COMMON_HEADERS)
	@echo "  Making dictionary for common headers/sources..."
	@rootcling -f $@ $^

#===============================================================
# clean-up and print command
#===============================================================

.PHONY: clean all

clean:
	@echo "  Cleaning..."
	@rm -f $(PACKAGES)
	@rm -f $(DICTIONARY)
	@rm -f $(DICTIONARY:%.cxx=%.h)
	@rm -f *rdict.pcm

.PHONY: print

print:
	@echo ${PACKAGES}
	@echo ${COMMON_HEADERS}
	@echo ${COMMON_SOURCES}
	@echo ${LDFLAGS}
	@echo ${CXXFLAGS}
