pkgdataexampledir = $(pkgdatadir)/example
dist_pkgdataexample_DATA = \
	$(srcdir)/*.aln \
	$(srcdir)/*.fasta* \
	$(srcdir)/*.pdb

dist-hook:
	rm -rf `find $(distdir) -name .svn`