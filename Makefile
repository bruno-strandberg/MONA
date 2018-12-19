include $(NMHDIR)/Makerules

#======================================================================================
# master make file that calls make in common_software/ fitter_software/ and for
# all directories in apps/
#======================================================================================

.PHONY: nmhlib fitlib apps

apps: nmhlib fitlib
	@for dir in apps/* ; do \
		$(MAKE) -C $${dir} ; \
	done

fitlib: nmhlib
	$(MAKE) -C $(FITSOFTDIR)

nmhlib:
	$(MAKE) -C $(NMHSOFTDIR)

.PHONY: clean

clean:
	$(MAKE) -C $(NMHSOFTDIR) clean
	$(MAKE) -C $(FITSOFTDIR) clean
	@for dir in apps/* ; do \
		$(MAKE) -C $${dir} clean ; \
	done
