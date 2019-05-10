include $(MONADIR)/Makerules

CFLAGS = $(CXXFLAGS) $(ROOTCFLAGS) $(JPPCFLAGS)
LFLAGS = $(ROOTLIBS) $(NMHLIB)

COMMON   = PIDAlpha.C PIDBeta.C PIDGamma.C PIDDelta.C
SOURCES  = $(filter-out $(COMMON), $(wildcard *.C))
PACKAGES = $(SOURCES:.C=)

all: $(PACKAGES)

$(PACKAGES): % : %.C $(COMMON)
	@echo "  Building $@ with pre-requisities $^..."
	@g++ -O3 -o $@ $^ $(LFLAGS) $(CFLAGS)


.PHONY: clean

clean:
	rm -f $(PACKAGES)


.PHONY: print

print:
	@echo $(COMMON)
	@echo $(SOURCES)
	@echo $(PACKAGES)
