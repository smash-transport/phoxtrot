/*
 * Main executable:
 * Compute the cross sections and print the results to separate output files.
 *
 */

#include "write_tables.cpp"

int main(int argc, char **argv) {

  std::string directory;

  if (argc == 1) {
    directory = '.';
  } else if (argc == 2) {
    directory = argv[1];
  } else {
    std::cout
        << "Too many arguments provided. Expected at most one argument: \n"
           "./phoxtrot PATH_TO_OUTPUT_DIRECTORY"
        << '\n';
    std::exit(0);
  }

  write_Total_Xsec_PiRho_PiGamma(directory);
  write_Diff_Xsec_PiRho_PiGamma(directory);
  write_Total_Xsec_PiPi_RhoGamma(directory);
  write_Diff_Xsec_PiPi_RhoGamma(directory);
}
