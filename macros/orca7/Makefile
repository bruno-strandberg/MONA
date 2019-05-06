
include $(MONADIR)/Makerules

CFLAGS = $(CXXFLAGS) $(ROOTCFLAGS) $(JPPCFLAGS)
LFLAGS = $(ROOTLIBS) $(ROOFITLIBS) $(NMHLIB) $(FITLIB)

# variable definition
#====================

# classes that are shared between applications; i.e. 
# all applications that are built will need the relevant sources
SHARED = ORCA7

SHEAD   = $(SHARED:=.h)
SSOURCE = $(SHARED:=.C)
DICT    = SharedDic.cxx
LDEF    = LinkDef.h

# sources of applications

PSOURCES  = contour.C
PROGRAMS  = $(PSOURCES:.C=)

# compilation
#=============

all: $(PROGRAMS)

$(PROGRAMS): % : %.C $(SSOURCE) $(DICT)
	@echo "  Building $@ with pre-requisities $^..."
	@g++ -O3 -o $@ $^ $(LFLAGS) $(CFLAGS)

$(DICT): $(SHEAD) $(LDEF)
	@echo "  Making dictionary $@..."
	@rootcling -f $@ $^

# clean-up and print command
#===========================

.PHONY: clean all

clean:
	@echo "  Cleaning..."
	@rm -f $(PROGRAMS)
	@rm -f $(DICT)
	@rm -f *rdict.pcm

.PHONY: print

print:
	@echo $(SSOURCE)
	@echo $(SHEAD)
	@echo $(DICT)
	@echo $(PSOURCES)
	@echo $(CFLAGS)
	@echo $(LFLAGS)
