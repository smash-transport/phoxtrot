/*
 * Main executable:
 * Compute the cross sections and print the results to a separate output.
 *
 */

#include "write_tables.cpp"

int main() {

  write_Total_Xsec_PiRho_PiGamma();
  write_Diff_Xsec_PiRho_PiGamma();
  write_Total_Xsec_PiPi_RhoGamma();
  write_Diff_Xsec_PiPi_RhoGamma();
}
