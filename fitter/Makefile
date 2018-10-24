
#===============================================================
# setup cxx and library flags
#===============================================================
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) 
ROOTLIBS     += -lRooFit -lRooFitCore

CXX           = g++
CXXFLAGS      = -g -Wall -std=c++11 -fPIC

NMHLIBS = -L $(NMHDIR)/common_software -lnmhsoft -L $(OSCPROBDIR) -lOscProb

CXXFLAGS  += $(ROOTCFLAGS)
LDFLAGS   += $(ROOTLIBS) $(NMHLIBS)

#===============================================================
# define variables for complilation
#===============================================================
PACKAGE1 = fitter
PACKAGE2 = asymmetry
HEADERS = $(filter-out Fitter.h, $(wildcard *.h))
SOURCES1 = $(filter-out Asymmetry.C, $(wildcard *.C))
SOURCES2 = $(filter-out Fitter.C, $(wildcard *.C))
TARGET1  = $(PACKAGE1)
TARGET2  = $(PACKAGE2)
DICTIONARY = dict.cxx

#for debugging
#debug:
#	@echo ${LDFLAGS}
#	@echo ${CXXFLAGS}
#	@echo ${SOURCES}
#	@echo ${HEADERS}

#===============================================================
# compilation
#===============================================================

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(DICTIONARY) $(SOURCES1)
	@echo "  Building $(PACKAGE1)..."
	@g++ -O3 -o $@ $^ $(LDFLAGS) $(CXXFLAGS)

$(TARGET2): $(DICTIONARY) $(SOURCES2)
	@echo "  Building $(PACKAGE2)..."
	@g++ -O3 -o $@ $^ $(LDFLAGS) $(CXXFLAGS)

$(DICTIONARY): $(HEADERS)
	@echo "  Making dictionary for $(PACKAGE)..."
	@rootcling -f $@ $^

clean:
	@echo "  Cleaning..."
	@rm -f $(TARGET1)
	@rm -f $(TARGET2)
	@rm -f $(DICTIONARY)
	@rm -f $(DICTIONARY:%.cxx=%.h)
	@rm -f *rdict.pcm

.PHONY: clean all

