
#=================================================================
# setup cxx flags and library flags
#=================================================================
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) 
ROOTLIBS     += -lRooFit -lRooFitCore

CXX           = g++
CXXFLAGS      = -g -Wall -std=c++11 -fPIC

CXXFLAGS  += $(ROOTCFLAGS)
LDFLAGS   += $(ROOTLIBS)

#=================================================================
# setup programs to be compiled
#=================================================================

PACKAGE1 = simfit
PACKAGE2 = compfit
HEADERS  = $(filter-out LinkDef.h, $(wildcard *.h)) #same header files for both 
SOURCES1 = $(filter-out Compfit.C, $(wildcard *.C))
SOURCES2 = $(filter-out Simfit.C, $(wildcard *.C))
TARGET1  = $(PACKAGE1)
TARGET2  = $(PACKAGE2)
DICTIONARY = dicfits.cxx
ROOTMAP    = libfits.rootmap

#for debugging
#debug:
#	@echo ${LDFLAGS}
#	@echo ${CXXFLAGS}
#	@echo ${SOURCES}
#	@echo ${HEADERS}

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(DICTIONARY) $(SOURCES1)
	@echo "  Building $(PACKAGE1)..."
	@g++ -O3 -o $@ $^ $(LDFLAGS) $(CXXFLAGS)

$(TARGET2): $(DICTIONARY) $(SOURCES2)
	@echo "  Building $(PACKAGE2)..."
	@g++ -O3 -o $@ $^ $(LDFLAGS) $(CXXFLAGS)

$(DICTIONARY): $(HEADERS) LinkDef.h
	@echo "  Making dictionary for..."
	@rootcling -f $@ -rml $(TARGET1) -rml $(TARGET2) -rmf $(ROOTMAP) $^

clean:
	@echo "  Cleaning..."
	@rm -f $(TARGET1)
	@rm -f $(TARGET2)
	@rm -f $(DICTIONARY)
	@rm -f $(DICTIONARY:%.cxx=%.h)
	@rm -f *rdict.pcm
	@rm -f $(ROOTMAP)

.PHONY: clean all

