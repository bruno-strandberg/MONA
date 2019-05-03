
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
DICT    = SharedDic
LDEF    = LinkDef.h

# sources of applications

PSOURCES  = contour.C
PROGRAMS  = $(PSOURCES:.C=)

# compilation
#=============

all: $(PROGRAMS)

# note that I get rid of the annoying dictionary file automatically; this is only required for compilation 
# to programs and can be removed, because no shared library object that requires the dict is created.
$(PROGRAMS): % : %.C $(SSOURCE) $(DICT).cxx
	@echo "  Building $@ with pre-requisities $^..."
	@g++ -O3 -o $@ $^ $(LFLAGS) $(CFLAGS) && rm $(DICT)*

$(DICT).cxx: $(SHEAD) $(LDEF)
	@echo "  Making dictionary $@..."
	@rootcling -f $@ $^

# clean-up and print command
#===========================

.PHONY: clean all

clean:
	@echo "  Cleaning..."
	@rm -f $(PROGRAMS)

.PHONY: print

print:
	@echo $(SSOURCE)
	@echo $(SHEAD)
	@echo $(DICT)
	@echo $(PSOURCES)
	@echo $(CFLAGS)
	@echo $(LFLAGS)
