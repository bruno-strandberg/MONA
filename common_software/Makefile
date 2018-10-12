#Compilation to ROOT library was a gift from Joao Coelho's OscProb package.

#get root CFLAGS and root LIBS
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) 

#set compiler and standard flags
CXX           = g++
CXXFLAGS      = -g -Wall -std=c++11 -fPIC

#add the ROOT flags, include dir is current dir
CXXFLAGS  += $(ROOTCFLAGS)
LDFLAGS   += $(ROOTLIBS)

# aanet use flag set in setenv.sh, headers are added to CPATH by setenv, no need for -I/...
ifeq ($(USE_AANET), true)
	CXXFLAGS += -DWAANET
        LDFLAGS  += -L${AADIR} -laa
endif

#for the rootcint command LinkDef.h needs to come last, that's why it is included
#only when Dictionary is created
PACKAGE    = nmhsoft
HEADERS    = $(filter-out LinkDef.h, $(wildcard *.h))
SOURCES    = $(wildcard *.C)
TARGET_LIB = lib$(PACKAGE).so
DICTIONARY = dic$(PACKAGE).cxx
ROOTMAP    = lib$(PACKAGE).rootmap

#for debugging
#debug:
#	@echo ${LDFLAGS}
#	@echo ${CXXFLAGS}
#	@echo ${SOURCES}
#	@echo ${HEADERS}

all: $(TARGET_LIB)

#target depends on dictionary (->depends on headers and linkdef) and sources
#$^ means prerequesities --> $^ = $HEADERS Linkdef.h $SOURCES
$(TARGET_LIB): $(DICTIONARY) $(SOURCES)
	@echo "  Building $(PACKAGE)..."
	@g++ -shared -O3 -o $@ $^ $(LDFLAGS) $(CXXFLAGS)

#$@ means file name of the target of the rule (dic$PACAGE)
#$^ means prerequesities --> HEADERS and Linkdef
$(DICTIONARY): $(HEADERS) LinkDef.h
	@echo "  Making dictionary for $(PACKAGE)..."
	@rootcling -f $@ -rml $(TARGET_LIB) -rmf $(ROOTMAP) $^

clean:
	@echo "  Cleaning $(PACKAGE)..."
	@rm -f $(TARGET_LIB)
	@rm -f $(DICTIONARY)
	@rm -f $(DICTIONARY:%.cxx=%.h)
	@rm -f *rdict.pcm
	@rm -f $(ROOTMAP)

.PHONY: clean all

