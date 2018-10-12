
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
PACKAGE = fitter
HEADERS = $(filter-out Fitter.h, $(wildcard *.h))
SOURCES = $(wildcard *.C)
TARGET  = $(PACKAGE)
DICTIONARY = dic$(PACKAGE).cxx

#for debugging
#debug:
#	@echo ${LDFLAGS}
#	@echo ${CXXFLAGS}
#	@echo ${SOURCES}
#	@echo ${HEADERS}

#===============================================================
# compilation
#===============================================================

all: $(TARGET)

$(TARGET): $(DICTIONARY) $(SOURCES)
	@echo "  Building $(PACKAGE)..."
	@g++ -O3 -o $@ $^ $(LDFLAGS) $(CXXFLAGS)

$(DICTIONARY): $(HEADERS)
	@echo "  Making dictionary for $(PACKAGE)..."
	@rootcling -f $@ $^

clean:
	@echo "  Cleaning $(PACKAGE)..."
	@rm -f $(TARGET)
	@rm -f $(DICTIONARY)
	@rm -f $(DICTIONARY:%.cxx=%.h)
	@rm -f *rdict.pcm

.PHONY: clean all

