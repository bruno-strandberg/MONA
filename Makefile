include $(MONADIR)/Makerules

#======================================================================================
# master make file that calls make in common_software/ fitter_software/ and for
# all directories in apps/
#======================================================================================

# default call of `make` builds the libraries and applications
#---------------------------------------------------------------------
.PHONY: all
all : apps

apps: nmhlib fitlib
	@for dir in apps/* ; do \
		$(MAKE) -C $${dir} || exit 1 ; \
	done

fitlib: nmhlib
	$(MAKE) -C $(FITSOFTDIR)

nmhlib: env
	$(MAKE) -C $(NMHSOFTDIR)

env:
ifeq ($(MONADIR),"")
	@echo "Environment variable MONADIR not set, please run 'source setenv.sh'. Exiting."
	exit 1
endif

# add a phony target to call the test-suite
#---------------------------------------------------------------------
.PHONY: test

test: tests
	python $(MONADIR)/run_tests.py

tests: all
	@for dir in $(shell ls -d tests/*/) ; do \
		$(MAKE) -C $${dir} || exit 1; \
	done

# add a phony target to call clean
#---------------------------------------------------------------------
.PHONY: clean

clean:
	$(MAKE) -C $(NMHSOFTDIR) clean
	$(MAKE) -C $(FITSOFTDIR) clean
	@for dir in apps/* ; do \
		$(MAKE) -C $${dir} clean ; \
	done
	@for dir in tests/* ; do \
		$(MAKE) -C $${dir} clean ; \
	done
