include $(MONADIR)/Makerules

CFLAGS = $(CXXFLAGS) $(ROOTCFLAGS) $(JPPCFLAGS)
LFLAGS = $(ROOTLIBS) $(NMHLIB)

SOURCES  = $(wildcard *.C)
HEADERS  = $(wildcard *.h)
PACKAGES = $(SOURCES:.C=)

all: $(PACKAGES)

$(PACKAGES): %: %.C
	@echo "  Building $@ with pre-requisities $^..."
	@g++ -O3 -o $@ $^ $(LFLAGS) $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(PACKAGES)

.PHONY: print

print:
	@echo $(SOURCES)
	@echo $(HEADERS)
	@echo $(PACKAGES)
	@echo $(LFLAGS)
	@echo $(CFLAGS)
