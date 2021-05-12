/*
 * Functions to print total and differential cross sections for all processes
 * to txt files.
 */

#include "include/common_functions.h"
#include "include/write_tables.h"
#include "include/TwoToTwo_differential_Xsections.h"
#include "include/TwoToTwo_total_Xsections.h"

void TwoToTwo_writing::write_Total_Xsec_PiRho_PiGamma(std::string dir) {
  std::ofstream Total_Xsections_PiRho_PiGamma;
  Total_Xsections_PiRho_PiGamma.open(dir + "/Total_Xsec_PiRho_PiGamma_mrho_" +
                                     std::to_string(int(mrho * 1000)) +
                                     "_MeV.txt");

  Total_Xsections_PiRho_PiGamma << "# sqrt(s), C11, C12, C13, C14, C15, C16"
                                << std::endl;

  double sigma_C11, sigma_C12, sigma_C13, sigma_C14, sigma_C15, sigma_C16;
  double sqrt_s;
  for (int i = 1; i < 1002; i += 1) {
    sqrt_s = mpion + mrho + i * 0.005;

    sigma_C11 = TwoToTwo_Tot_Xsections::total_xsection_C11(sqrt_s * sqrt_s);
    sigma_C12 = TwoToTwo_Tot_Xsections::total_xsection_C12(sqrt_s * sqrt_s);
    sigma_C13 = TwoToTwo_Tot_Xsections::total_xsection_C13(sqrt_s * sqrt_s);
    sigma_C14 = TwoToTwo_Tot_Xsections::total_xsection_C14(sqrt_s * sqrt_s);
    sigma_C15 = TwoToTwo_Tot_Xsections::total_xsection_C15(sqrt_s * sqrt_s);
    sigma_C16 = TwoToTwo_Tot_Xsections::total_xsection_C16(sqrt_s * sqrt_s);

    Total_Xsections_PiRho_PiGamma
        << sqrt_s << "\t" << sigma_C11 << "\t" << sigma_C12 << "\t" << sigma_C13
        << "\t" << sigma_C14 << "\t" << sigma_C15 << "\t" << sigma_C16
        << std::endl;
  }
  Total_Xsections_PiRho_PiGamma.close();
}

void TwoToTwo_writing::write_Diff_Xsec_PiRho_PiGamma(std::string dir) {
  std::ofstream Diff_Xsections_PiRho_PiGamma;
  Diff_Xsections_PiRho_PiGamma.open(dir + "/Diff_Xsec_PiRho_PiGamma_mrho_" +
                                    std::to_string(int(mrho * 1000)) +
                                    "_MeV.txt");

  double sqrt_s = 1.0;
  Diff_Xsections_PiRho_PiGamma
      << "# Differential sigma at sqrt(s) = " + std::to_string(sqrt_s) +
             " GeV as a function of mandelstam t and the scattering angle "
             "\u03F4"
      << std::endl;
  Diff_Xsections_PiRho_PiGamma
      << "# t, \u03F4 (degree), C11, C12, C13, C14, C15, C16" << std::endl;

  double tmin = min_mandelstam_t(sqrt_s * sqrt_s, mpion, mrho, mpion);
  double tmax = max_mandelstam_t(sqrt_s * sqrt_s, mpion, mrho, mpion);
  double increment = (tmax - tmin) / 1000.0;
  double t, theta;
  double diff_sigma_C11, diff_sigma_C12, diff_sigma_C13;
  double diff_sigma_C14, diff_sigma_C15, diff_sigma_C16;

  for (int i = 0; i < 1001; i += 1) {
    t = tmin + i * increment;
    theta = theta_from_t(t, sqrt_s, mpion, mrho, mpion);

    diff_sigma_C11 = TwoToTwo_Diff_Xsections::diff_xsection_C11(t, sqrt_s * sqrt_s);
    diff_sigma_C12 = TwoToTwo_Diff_Xsections::diff_xsection_C12(t, sqrt_s * sqrt_s);
    diff_sigma_C13 = TwoToTwo_Diff_Xsections::diff_xsection_C13(t, sqrt_s * sqrt_s);
    diff_sigma_C14 = TwoToTwo_Diff_Xsections::diff_xsection_C14(t, sqrt_s * sqrt_s);
    diff_sigma_C15 = TwoToTwo_Diff_Xsections::diff_xsection_C15(t, sqrt_s * sqrt_s);
    diff_sigma_C16 = TwoToTwo_Diff_Xsections::diff_xsection_C16(t, sqrt_s * sqrt_s);

    Diff_Xsections_PiRho_PiGamma
        << t << "\t" << theta << "\t" << diff_sigma_C11 << "\t"
        << diff_sigma_C12 << "\t" << diff_sigma_C13 << "\t" << diff_sigma_C14
        << "\t" << diff_sigma_C15 << "\t" << diff_sigma_C16 << std::endl;
  }
  Diff_Xsections_PiRho_PiGamma.close();
}

void TwoToTwo_writing::write_Total_Xsec_PiPi_RhoGamma(std::string dir) {
  std::ofstream Total_Xsections_PiPi_RhoGamma;
  Total_Xsections_PiPi_RhoGamma.open(dir + "/Total_Xsec_PiPi_RhoGamma_mrho_" +
                                     std::to_string(int(mrho * 1000)) +
                                     "_MeV.txt");
  Total_Xsections_PiPi_RhoGamma << "# sqrt(s), C21, C22" << std::endl;

  double sigma_C21, sigma_C22;
  double sqrt_s;
  for (int i = 1; i < 1002; i += 1) {
    // Threshold depends on mass of rho meson
    if (mrho > 2 * mpion) {
      sqrt_s = mrho + i * 0.005;
    } else {
      sqrt_s = mpion + mpion + i * 0.005;
    }

    sigma_C21 = TwoToTwo_Tot_Xsections::total_xsection_C21(sqrt_s * sqrt_s);
    sigma_C22 = TwoToTwo_Tot_Xsections::total_xsection_C22(sqrt_s * sqrt_s);

    Total_Xsections_PiPi_RhoGamma << sqrt_s << "\t" << sigma_C21 << "\t"
                                  << sigma_C22 << std::endl;
  }
  Total_Xsections_PiPi_RhoGamma.close();
}

void TwoToTwo_writing::write_Diff_Xsec_PiPi_RhoGamma(std::string dir) {
  std::ofstream Diff_Xsections_PiPi_RhoGamma;
  Diff_Xsections_PiPi_RhoGamma.open(dir + "/Diff_Xsec_PiPi_RhoGamma_mrho_" +
                                    std::to_string(int(mrho * 1000)) +
                                    "_MeV.txt");

  double sqrt_s = 1.0;
  Diff_Xsections_PiPi_RhoGamma
      << "# Differential sigma at sqrt(s) = " + std::to_string(sqrt_s) +
             " GeV as a function of mandelstam t and the scattering angle "
             "\u03F4"
      << std::endl;
  Diff_Xsections_PiPi_RhoGamma << "# t, \u03F4 (degree), C21, C22" << std::endl;

  double tmin = min_mandelstam_t(sqrt_s * sqrt_s, mpion, mpion, mrho);
  double tmax = max_mandelstam_t(sqrt_s * sqrt_s, mpion, mpion, mrho);
  double increment = (tmax - tmin) / 1000.0;
  double t, theta;
  double diff_sigma_C21, diff_sigma_C22;

  for (int i = 0; i < 1001; i += 1) {
    t = tmin + i * increment;
    theta = theta_from_t(t, sqrt_s, mpion, mpion, mrho);

    diff_sigma_C21 = TwoToTwo_Diff_Xsections::diff_xsection_C21(t, sqrt_s * sqrt_s);
    diff_sigma_C22 = TwoToTwo_Diff_Xsections::diff_xsection_C22(t, sqrt_s * sqrt_s);

    Diff_Xsections_PiPi_RhoGamma << t << "\t" << theta << "\t" << diff_sigma_C21
                                 << "\t" << diff_sigma_C22 << "\t" << std::endl;
  }
  Diff_Xsections_PiPi_RhoGamma.close();
}
