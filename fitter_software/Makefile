
include $(NMHDIR)/Makerules

#===============================================================
# setup cxx and library flags from global variables
#===============================================================

CFLAGS = $(CXXFLAGS) $(ROOTCFLAGS)
LFLAGS = $(ROOTLIBS) $(ROOFITLIBS) $(NMHLIB) $(OSCLIB)

#===============================================================
# define variables for complilation
#===============================================================
HEADERS = $(wildcard *.h)
SOURCES = $(wildcard *.C)
OBJ     = lib$(FITSOFT).so
DICT    = $(FITSOFT).cxx

#===============================================================
# compilation
#===============================================================

all: $(OBJ)

$(OBJ): $(SOURCES) $(DICT)
	@echo "  Building $@ with pre-requisities $^..."
	@g++ -shared -O3 -o $@ $^ $(LFLAGS) $(CFLAGS)

$(DICT): $(HEADERS)
	@echo "  Making dictionary $@..."
	@rootcling -f $@ $^

#===============================================================
# clean-up and print command
#===============================================================

.PHONY: clean all

clean:
	@echo "  Cleaning..."
	@rm -f $(OBJ)
	@rm -f $(DICT)
	@rm -f *rdict.pcm

.PHONY: print

print:
	@echo ${HEADERS}
	@echo ${SOURCES}
	@echo ${OBJ}
	@echo ${CFLAGS}
