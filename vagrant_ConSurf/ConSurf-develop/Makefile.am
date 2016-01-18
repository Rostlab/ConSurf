dist_bin_SCRIPTS = consurf
dist_pkgdata_DATA = \
	CONSURF_FUNCTIONS.pm \
	CONSURF_CONSTANTS.pm \
	MSA_parser.pm \
	pdbParser.pm \
	prepareMSA.pm \
	TREE_parser.pm  \
	cp_rasmol_gradesPE_and_pipe.pm \
	consurfrc.default \
	consurf_new.py \
	chimera_consurf.cmd 


dist_pkgdata_SCRIPTS =

# lkajan: CentOS 5 autoconf does not define docdir
docdir = $(datadir)/doc/$(PACKAGE)
dist_doc_DATA = AUTHORS README

SUBDIRS = example output

man_MANS = consurf.1

consurf.1: consurf
	sed -e 's|__datadir__|$(datadir)|g;s|__docdir__|$(docdir)|g;s|__pkgdatadir__|$(pkgdatadir)|g;s|__PREFIX__|$(prefix)|g;s|__sysconfdir__|$(sysconfdir)|g;s|__PACKAGE_VERSION__|$(VERSION)|g;' "$<" | \
	pod2man -c 'User Commands' -r "$(VERSION)" -name $(shell tr '[:lower:]' '[:upper:]' <<< "$(basename $@)") > "$@"

clean-local:
	rm -f $(man_MANS)

install-data-hook:
	for f in '$(DESTDIR)$(pkgdatadir)/consurfrc.default'; do \
		sed -i -e 's|__datadir__|$(datadir)|g;s|__docdir__|$(docdir)|g;s|__pkgdatadir__|$(pkgdatadir)|g;s|__PREFIX__|$(prefix)|g;s|__sysconfdir__|$(sysconfdir)|g;s|__PACKAGE_VERSION__|$(VERSION)|g;' "$$f"; \
	done

install-exec-hook:
	sed -i -e 's|__datadir__|$(datadir)|g;s|__docdir__|$(docdir)|g;s|__pkgdatadir__|$(pkgdatadir)|g;s|__PREFIX__|$(prefix)|g;s|__sysconfdir__|$(sysconfdir)|g;s|__PACKAGE_VERSION__|$(VERSION)|g;' "$(DESTDIR)$(bindir)/consurf"

