include $(NMHDIR)/Makerules

#======================================================================================
# master make file that calls make in common_software/ fitter_software/ and for
# all directories in apps/
#======================================================================================

# default call of `make` builds the libraries, applications and tests
#---------------------------------------------------------------------
.PHONY all: tests

tests: apps
	@for dir in tests/* ; do \
		$(MAKE) -C $${dir} ; \
	done

apps: nmhlib fitlib
	@for dir in apps/* ; do \
		$(MAKE) -C $${dir} ; \
	done

fitlib: nmhlib
	$(MAKE) -C $(FITSOFTDIR)

nmhlib:
	$(MAKE) -C $(NMHSOFTDIR)

# add a phony target to call the test-suite
#---------------------------------------------------------------------
.PHONY: test

test: tests
	python $(NMHDIR)/run_tests.py

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
