include $(MONADIR)/Makerules

CFLAGS = $(CXXFLAGS) $(ROOTCFLAGS) $(JPPCFLAGS) -I./parsers/
LFLAGS = $(ROOTLIBS) $(NMHLIB)

PARSERS  = $(wildcard parsers/*.C)
SOURCES  = $(wildcard *.C)
APPS     = $(SOURCES:.C=)

all: $(APPS)

$(APPS): % : %.C $(PARSERS)
	@echo "  Building $@..."
	@g++ -O3 -o $@ $^ $(LFLAGS) $(CFLAGS)


.PHONY: clean

clean:
	rm -f $(APPS)


.PHONY: print

print:
	@echo $(PARSERS)
	@echo $(SOURCES)
	@echo $(APPS)
	@echo $(CFLAGS)
	@echo $(LFLAGS)
