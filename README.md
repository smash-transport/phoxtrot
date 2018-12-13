# PHOXTROT

**Caution: This repository is still very much work in progress...**

Welcome to PHOXTROT (**PHO**ton **X**sections for **TR**ansp**O**r**T**).
This repository contains a collection of cross sections useful to describe the
production of photons in a hadronic medium. See * **paper(once public)** * for
further details.

The presented cross sections and processes rely on a theoretical framework
described by S. Turbide, R. Rapp and C. Gale in [
Phys. Rev. C 69, 014903 (2004)](https://arxiv.org/pdf/hep-ph/0308085.pdf).
The matrix elements are taken from the
[PhD thesis](http://digitool.library.mcgill.ca/R/?func=dbin-jump-full&object_id=102221&local_base=GEN01-MCG02)
of Simon Turbide.

In this repository, the differential and total cross sections are presented in
C++ readable format. The entire repository is wrapped by a CMake routine to
easily produce tables and plots for the given processes. To build PHOXTROT and
create the executable, use the following commands:

    mkdir build
    cd build
    cmake ..
    make phoxtrot

The executable `phoxtrot` is then located in the build directory. Upon execution, txt files for the differential and total cross section are created that contain the computed cross sections for varying values of t and s. To execute phoxtrot and to create the cross section files, type

    make phoxtrot_tables

The generated cross section files can be found in `PATH_TO_BUILD/tables`.
