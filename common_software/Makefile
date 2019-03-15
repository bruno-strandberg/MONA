
include $(MONADIR)/Makerules

#===============================================================
# setup cxx and library flags from global variables
#===============================================================

CFLAGS = $(CXXFLAGS) $(ROOTCFLAGS) $(AACFLAGS)
LFLAGS = $(ROOTLIBS) $(AALIB)

#===============================================================
# define variables for complilation
#===============================================================

HEADERS = $(filter-out LinkDef.h, $(wildcard *.h))
SOURCES = $(wildcard *.C)
OBJ     = lib$(NMHSOFT).so
DICT    = dic$(NMHSOFT).cxx
ROOTMAP = dic$(NMHSOFT).rootmap

#===============================================================
# compilation
#===============================================================

all: $(OBJ)

$(OBJ): $(SOURCES) $(DICT)
	@echo "  Building $@ with pre-requisities $^..."
	@g++ -shared -O3 -o $@ $^ $(LFLAGS) $(CFLAGS)

$(DICT): $(HEADERS) LinkDef.h
	@echo "  Making dictionary $@..."
	@rootcling -f $@ -rml $(OBJ) -rmf $(ROOTMAP) $^

#===============================================================
# clean-up and print command
#===============================================================

.PHONY: clean all

clean:
	@echo "  Cleaning..."
	@rm -f $(OBJ)
	@rm -f $(DICT)
	@rm -f $(ROOTMAP)
	@rm -f *rdict.pcm

.PHONY: print

print:
	@echo ${HEADERS}
	@echo ${SOURCES}
	@echo ${OBJ}
	@echo ${DICT}
	@echo ${ROOTMAP}
	@echo ${CFLAGS}
	@echo ${LFLAGS}
