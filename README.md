# PHOXTROT

Welcome to PHOXTROT (**PHO**ton **X**sections for **TR**ansp**O**r**T**).
This repository contains a collection of cross sections useful to describe the
production of photons in a hadronic medium. See for example [Phys. Rev. D 99, 114021 (2019)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.99.114021) for
further details.

PHOXTROT contains cross sections for two types of processes:
- **2 <-> 2 Scatterings**: π + π -> ρ + γ and π + ρ -> π + γ
- **Bremsstrahlung**: π + π -> π + π + γ

The presented cross sections for 2 <-> 2 scatterings rely on a theoretical framework
described by S. Turbide, R. Rapp and C. Gale in [
Phys. Rev. C 69, 014903 (2004)](https://arxiv.org/pdf/hep-ph/0308085.pdf).
The matrix elements are taken from the
[PhD thesis](http://digitool.library.mcgill.ca/R/?func=dbin-jump-full&object_id=102221&local_base=GEN01-MCG02)
of Simon Turbide. <br>
The presented cross sections for bremsstrahlung processes are calculated within the one-boson exchange model.

In this repository, the differential and total cross sections are presented in
C++ readable format. The entire repository is wrapped by a CMake routine to
easily produce tables and plots for the given processes.

##### Prerequisites
To properly use PHOXTROT, make sure the following prerequisites are met:
- cmake >= 3.15.4
- the GNU Scientific Library >= 2.0
- Python modules: matplotlib, numpy, seaborn, argparse, os, sys, glob

##### Compilation
To build PHOXTROT and create the executable, use the following commands:

    mkdir build
    cd build
    cmake ..
    make phoxtrot

The executable `phoxtrot` is located in the build directory. Upon
execution, txt files for the differential and total cross section are created
that contain the computed cross sections for varying values of theta, k and s. To
execute phoxtrot, create the cross section files and plot them, type

    make plots

The generated cross section files can then be found in `PATH_TO_BUILD/tables`
and the corresponding plots are located in `PATH_TO_BUILD/plots`. Note that the tables and plots for all differential cross sections are created at a fixed value of √s = 1.0 GeV. If necessary, this value may be changed in the functions `TwoToTwo_writing::write_Diff_Xsec_PiRho_PiGamma`, `TwoToTwo_writing::write_Diff_Xsec_PiPi_RhoGamma`, `Brems_writing::write_dSigma_dk`, and `Brems_writing::write_dSigma_dTheta`.
