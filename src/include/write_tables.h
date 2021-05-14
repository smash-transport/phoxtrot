// write_tables.h

#ifndef WRITE_TABLES_H
#define WRITE_TABLES_H

#include <string>

class TwoToTwo_writing {

  public:
    // Write cross sections for (π + ρ0 ->  π + γ) processes
    static void write_Total_Xsec_PiRho_PiGamma(std::string dir);
    static void write_Diff_Xsec_PiRho_PiGamma(std::string dir);

    // Write cross sections for (π + π -> ρ0 + γ) processes
    static void write_Total_Xsec_PiPi_RhoGamma(std::string dir);
    static void write_Diff_Xsec_PiPi_RhoGamma(std::string dir);
};

class Brems_writing {

  public:
    // Write total cross sections
    static void write_Total_Xsec(std::string dir);

    // Write differential cross sections for dSigma/dk
    static void write_dSigma_dk(std::string dir);

    // Write differential cross sections for dSigma/dtheta
    static void write_dSigma_dTheta(std::string dir);
};

#endif
