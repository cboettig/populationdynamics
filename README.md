[![DOI](https://zenodo.org/badge/3437830.svg)](https://zenodo.org/badge/latestdoi/3437830)



- Package: populationdynamics
- Title: An R package with tools to simulate various population dynamics models in ecology
- Version: 0.0-1
- Date: 2012-02-13
- Author: Carl Boettiger <cboettig@gmail.com>
- Maintainer: Carl Boettiger <cboettig@gmail.com>
- Description: Simulation tools, largely for individual based models
   using the Gillespie algorithm and fast C code.
- License: BSD
- URL: https://github.com/cboettig/populationdynamics
- BugReports: https://github.com/cboettig/populationdynamics/issues
- SystemRequirements: Gnu Scientific Library version >= 1.8

---

## Notes 

This package has been developed to house functions that I commonly use
in my own research.  It is not intended as a general-purpose toolbox
and thus may not have the polish or robustness you might expect from an
R package intended for the use of others, such as you find on CRAN.

These functions are provided as an R package only to provide a module for
my own various research projects. Nonetheless, the package is licensed
under permissive BSD license and you are more then welcome to adapt it as
you will.  The most successful reuse may involve extracting particular
bits of code, such as the Gillespie C engine, to extend rather than
installing the package whole-cloth.

Many functions in this toolbox make use of fast, parallelized
(OpenMP) C code for improved performance of the exact, continuous time
individual-based stochastic simulations using the Gillespie algorithm,
including the generation of many replicate simulations.  Unfortunately,
this adds additional hurdles to installation and portability of the
package.

Note that the package requires the GNU Scientific C Library be
installed to run.  The package will attempt to find and link the libary
automatically using autoconf tools, though sometimes manual intervention
is required. See R manuals on installing R packages with autoconf for
more details.  The package has been tested exclusively on Linux systems,
including the high performance computing facility Carver operated by
the NERSC supercomputing center in Berkeley.


## Installation hints

One Mac user found it necessary to add a symbolic link to the appropriate C compiler to successfully install this package on Maverick. See https://github.com/cboettig/prosecutors-fallacy/issues/2#issue-27830362
