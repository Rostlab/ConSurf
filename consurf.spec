Summary: Identification of Functional Regions in Proteins
Name: consurf
Version: 1.0.0
Release: 1
License: NON COMMERCIAL SOFTWARE LICENSE AGREEMENT
Group: Applications/Science
Source: ftp://rostlab.org/%{name}/%{name}-%{version}.tar.gz
URL: http://consurf.tau.ac.il/
BuildArch: noarch
BuildRoot: %{_tmppath}/%{name}-%{version}-root
BuildRequires: perl
Requires: blast2, clustalw, muscle, rate4site

%description
 The ConSurf system identifies functional regions in proteins. It is mainly based on the Rate4Site algorithm
 for detecting conserved amino-acid sites by computing the relative evolutionary rate for
 each site in the multiple sequence alignment (MSA).
 .
 The ConSurf work in several modes
 .
 1. Given Protein PDB File.
 2. Given Multiple Sequence Alignment (MSA) and Protein PDB File.
 3. Given Multiple Sequence Alignment (MSA), Phylogenetic Tree, and PDB File.
 .
 The script is using the user provided MSA (and Phylogenetic tree if available) to calculate the conservation score for each position in the MSA based on the Rate4Site algorithm (Mayrose, I.,
 Graur, D., Ben-Tal, N., and Pupko, T. 2004. Comparison of site-specific rate-inference methods: Bayesian methods are superior. Mol Biol Evol 21: 1781-1791).
 .
 When running in the first mode the scripts the MSA is automatically build the MSA for the given protein based on ConSurf protocol.

%prep
%setup -q

%global __perl_requires %{_builddir}/%{name}-%{version}/%{name}-req
chmod +x %{__perl_requires}

%build
%configure
make

%install
rm -rf $RPM_BUILD_ROOT
make DESTDIR=${RPM_BUILD_ROOT} install

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%doc %{_docdir}/%{name}/AUTHORS
%doc %{_docdir}/%{name}/README
%doc COPYING
%{_bindir}/*
%{_mandir}/*/*
%{_datadir}/%{name}/*

%changelog
* Tue Sep 27 2011 Guy Yachdav <gyachdav@rostlab.org> - 1.0.0-1
- new upstream
- First rpm package
