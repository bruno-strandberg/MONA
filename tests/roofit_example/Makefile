
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
HEADERS  = $(wildcard *.h)                             #same header files for both (Compfit/Simfit only .C file)
SOURCES1 = $(filter-out Compfit.C, $(wildcard *.C))
SOURCES2 = $(filter-out Simfit.C, $(wildcard *.C))
TARGET1  = $(PACKAGE1)
TARGET2  = $(PACKAGE2)
DICTIONARY = dicfits.cxx

#for debugging
#debug:
#	@echo ${LDFLAGS}
#	@echo ${CXXFLAGS}
#	@echo ${SOURCES}
#	@echo ${HEADERS}

#=================================================================
# build commands
#=================================================================

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(DICTIONARY) $(SOURCES1)
	@echo "  Building $(PACKAGE1)..."
	@g++ -O3 -o $@ $^ $(LDFLAGS) $(CXXFLAGS)

$(TARGET2): $(DICTIONARY) $(SOURCES2)
	@echo "  Building $(PACKAGE2)..."
	@g++ -O3 -o $@ $^ $(LDFLAGS) $(CXXFLAGS)

$(DICTIONARY): $(HEADERS)
	@echo "  Making the dictionary..."
	@rootcling -f $@ $^

clean:
	@echo "  Cleaning..."
	@rm -f $(TARGET1)
	@rm -f $(TARGET2)
	@rm -f $(DICTIONARY)
	@rm -f $(DICTIONARY:%.cxx=%.h)
	@rm -f *rdict.pcm

.PHONY: clean all

