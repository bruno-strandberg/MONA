
#get root CFLAGS and root LIBS
ROOTCFLAGS    = $(shell root-config --cflags) -fPIC
ROOTLIBS      = $(shell root-config --nonew --libs) 

#set compiler and standard flags
CXX           = g++
CXXFLAGS      = -g -Wall -std=c++11

#add the ROOT flags, include dir is current dir
CXXFLAGS  += $(ROOTCFLAGS)
LDFLAGS   += $(ROOTLIBS)

#for the rootcint command LinkDef.h needs to come last, that's why it is included
#only when Dictionary is created
PACKAGE = commonsoft
HEADERS = $(filter-out $(CURDIR)/LinkDef.h, $(wildcard $(CURDIR)/*.h))
SOURCES = $(wildcard $(CURDIR)/*.C)
TARGET_LIB = $(CURDIR)/lib$(PACKAGE).so
DICTIONARY = $(CURDIR)/dic$(PACKAGE).cxx

all: $(TARGET_LIB)

#target depends on dictionary (->depends on headers and linkdef) and sources
#$^ means prerequesities --> $^ = $HEADERS Linkdef.h $SOURCES
$(TARGET_LIB): $(DICTIONARY) $(SOURCES)
	@echo "  Building $(PACKAGE)..."
	@g++ -shared -O3 -o$@ $(LDFLAGS) $(CXXFLAGS) $^

#$@ means file name of the target of the rule (dic$PACAGE)
#$^ means prerequesities --> HEADERS and Linkdef
$(DICTIONARY): $(HEADERS) LinkDef.h
	@echo "  Making dictionary for $(PACKAGE)..."
	@rootcint -f $@ -c $^

clean:
	@echo "  Cleaning $(PACKAGE)..."
	@rm -f $(TARGET_LIB)
	@rm -f $(DICTIONARY)
	@rm -f $(DICTIONARY:%.cxx=%.h)

.PHONY: clean all

