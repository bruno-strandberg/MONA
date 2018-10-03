#Compilation to ROOT library was a gift from Joao Coelho's OscProb package.

#get root CFLAGS and root LIBS
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) 
ROOTLIBS     += -lRooFit -lRooFitCore

#set compiler and standard flags
CXX           = g++
CXXFLAGS      = -g -Wall -std=c++11 -fPIC

#add the ROOT flags, include dir is current dir
CXXFLAGS  += $(ROOTCFLAGS)
LDFLAGS   += $(ROOTLIBS)

#for the rootcint command LinkDef.h needs to come last, that's why it is included
#only when Dictionary is created
PACKAGE = fit
HEADERS = $(filter-out $(CURDIR)/LinkDef.h, $(wildcard $(CURDIR)/*.h))
SOURCES = $(wildcard $(CURDIR)/*.C)
TARGET  = $(CURDIR)/$(PACKAGE)
DICTIONARY = $(CURDIR)/dic$(PACKAGE).cxx
ROOTMAP    = $(CURDIR)/lib$(PACKAGE).rootmap

#for debugging
#debug:
#	@echo ${LDFLAGS}
#	@echo ${CXXFLAGS}
#	@echo ${SOURCES}
#	@echo ${HEADERS}

all: $(TARGET)

#target depends on dictionary (->depends on headers and linkdef) and sources
#$^ means prerequesities --> $^ = $HEADERS Linkdef.h $SOURCES
$(TARGET): $(DICTIONARY) $(SOURCES)
	@echo "  Building $(PACKAGE)..."
	@g++ -O3 -o $@ $^ $(LDFLAGS) $(CXXFLAGS)

#$@ means file name of the target of the rule (dic$PACAGE)
#$^ means prerequesities --> HEADERS and Linkdef
$(DICTIONARY): $(HEADERS) LinkDef.h
	@echo "  Making dictionary for $(PACKAGE)..."
	@rootcling -f $@ -rml $(TARGET) -rmf lib$(PACKAGE).rootmap $^

clean:
	@echo "  Cleaning $(PACKAGE)..."
	@rm -f $(TARGET)
	@rm -f $(DICTIONARY)
	@rm -f $(DICTIONARY:%.cxx=%.h)
	@rm -f *rdict.pcm
	@rm -f $(ROOTMAP)

.PHONY: clean all

