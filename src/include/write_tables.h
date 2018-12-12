// write_tables.h

#ifndef WRITE_TABLES_H
#define WRITE_TABLES_H

// Write cross sections for (π + ρ0 ->  π + γ) processes
void write_Total_Xsec_PiRho_PiGamma();
void write_Diff_Xsec_PiRho_PiGamma();

// Write cross sections for (π + π -> ρ0 + γ) processes
void write_Total_Xsec_PiPi_RhoGamma();
void write_Diff_Xsec_PiPi_RhoGamma();

// Declarations of variables

double sigma_C11;
double sigma_C12;
double sigma_C13;
double sigma_C14;
double sigma_C15;
double sigma_C16;
double sigma_C21;
double sigma_C22;

double diff_sigma_C11;
double diff_sigma_C12;
double diff_sigma_C13;
double diff_sigma_C14;
double diff_sigma_C15;
double diff_sigma_C16;
double diff_sigma_C21;
double diff_sigma_C22;

#endif
