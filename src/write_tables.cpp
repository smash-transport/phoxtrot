/*
 * Functions to print total and differential cross sections for all processes
 * to txt files.
 */

#include<vector>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline2d.h>

#include "include/common_functions.h"
#include "include/write_tables.h"
#include "include/TwoToTwo_differential_Xsections.h"
#include "include/TwoToTwo_total_Xsections.h"
#include "include/bremsstrahlung_Xsections.h"

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

void Brems_writing::write_Total_Xsec(std::string dir) {

  // Read in tabularized total cross sections as vectors
  std::vector<double> sqrts_vec = BREMS_SQRTS;
  std::vector<double> sigma_pipi_pipi_opp_vec = BREMS_PIPI_PIPI_OPP_SIG;
  std::vector<double> sigma_pipi_pipi_same_vec = BREMS_PIPI_PIPI_SAME_SIG;
  std::vector<double> sigma_pipi0_pipi0_vec = BREMS_PIPI0_PIPI0_SIG;
  std::vector<double> sigma_pipi_pi0pi0_vec = BREMS_PIPI_PI0PI0_SIG;
  std::vector<double> sigma_pi0pi0_pipi_vec = BREMS_PI0PI0_PIPI_SIG;

  // Cast vectors into arrays such that gsl can cope with them
  const double* sqrts = &sqrts_vec[0];
  const double* sigma_pipi_pipi_opp = &sigma_pipi_pipi_opp_vec[0];
  const double* sigma_pipi_pipi_same = &sigma_pipi_pipi_same_vec[0];
  const double* sigma_pipi0_pipi0 = &sigma_pipi0_pipi0_vec[0];
  const double* sigma_pipi_pi0pi0 = &sigma_pipi_pi0pi0_vec[0];
  const double* sigma_pi0pi0_pipi = &sigma_pi0pi0_pipi_vec[0];

  // Get number of entries
  int length = sqrts_vec.size();

  // Set interpolations:
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline_pipi_pipi_opp = gsl_spline_alloc(gsl_interp_linear, length);
  gsl_spline_init (spline_pipi_pipi_opp, sqrts, sigma_pipi_pipi_opp, length);

  gsl_spline *spline_pipi_pipi_same = gsl_spline_alloc(gsl_interp_linear, length);
  gsl_spline_init (spline_pipi_pipi_same, sqrts, sigma_pipi_pipi_same, length);

  gsl_spline *spline_pipi0_pipi0 = gsl_spline_alloc(gsl_interp_linear, length);
  gsl_spline_init (spline_pipi0_pipi0, sqrts, sigma_pipi0_pipi0, length);

  gsl_spline *spline_pipi_pi0pi0 = gsl_spline_alloc(gsl_interp_linear, length);
  gsl_spline_init (spline_pipi_pi0pi0, sqrts, sigma_pipi_pi0pi0, length);

  gsl_spline *spline_pi0pi0_pipi = gsl_spline_alloc(gsl_interp_linear, length);
  gsl_spline_init (spline_pi0pi0_pipi, sqrts, sigma_pi0pi0_pipi, length);


  double sqrts_val;
  double sigma_val_pipi_pipi_opp, sigma_val_pipi_pipi_same, sigma_val_pipi0_pipi0, sigma_val_pipi_pi0pi0, sigma_val_pi0pi0_pipi;

  std::ofstream Brems_Total_Xsections;
  Brems_Total_Xsections.open(dir + "/Brems_Total_Xsections.txt");

  Brems_Total_Xsections << "# sqrt(s), pi+- pi-+ -> pi+- pi-+ gamma, pi+- pi+- -> pi+- pi+- gamma, pi pi0 -> pi pi0 gamma, pi pi -> pi0 pi0 gamma, pi0 pi0 -> pi pi gamma"
                                << std::endl;
  Brems_Total_Xsections << "# sqrt(s) in GeV, Sigma in mbarn" << std::endl;


  for (int i = 0; i < 470; i += 1) {
    sqrts_val = 0.3 + 0.01 * i;
    sigma_val_pipi_pipi_opp = gsl_spline_eval(spline_pipi_pipi_opp, sqrts_val, acc);
    sigma_val_pipi_pipi_same = gsl_spline_eval(spline_pipi_pipi_same, sqrts_val, acc);
    sigma_val_pipi0_pipi0 = gsl_spline_eval(spline_pipi0_pipi0, sqrts_val, acc);
    sigma_val_pipi_pi0pi0 = gsl_spline_eval(spline_pipi_pi0pi0, sqrts_val, acc);
    sigma_val_pi0pi0_pipi = gsl_spline_eval(spline_pi0pi0_pipi, sqrts_val, acc);

    Brems_Total_Xsections << sqrts_val << "\t" << sigma_val_pipi_pipi_opp
        << "\t" << sigma_val_pipi_pipi_same << "\t" << sigma_val_pipi0_pipi0
        << "\t" << sigma_val_pipi_pi0pi0 << "\t" << sigma_val_pi0pi0_pipi << std::endl;
  }
  Brems_Total_Xsections.close();

  gsl_spline_free (spline_pipi_pipi_opp);
  gsl_spline_free (spline_pipi_pipi_same);
  gsl_spline_free (spline_pipi0_pipi0);
  gsl_spline_free (spline_pipi_pi0pi0);
  gsl_spline_free (spline_pi0pi0_pipi);
  gsl_interp_accel_free (acc);
}

void Brems_writing::write_dSigma_dk(std::string dir) {
  // Read in tabularized differential cross sections as vectors
  std::vector<double> sqrts_vec = BREMS_SQRTS;
  std::vector<double> photon_momentum_vec = BREMS_K;
  std::vector<double> dsigma_dk_pipi_pipi_opp_vec = BREMS_PIPI_PIPI_OPP_DIFF_SIG_K;
  std::vector<double> dsigma_dk_pipi_pipi_same_vec = BREMS_PIPI_PIPI_SAME_DIFF_SIG_K;
  std::vector<double> dsigma_dk_pipi0_pipi0_vec = BREMS_PIPI0_PIPI0_DIFF_SIG_K;
  std::vector<double> dsigma_dk_pipi_pi0pi0_vec = BREMS_PIPI_PI0PI0_DIFF_SIG_K;
  std::vector<double> dsigma_dk_pi0pi0_pipi_vec = BREMS_PI0PI0_PIPI_DIFF_SIG_K;

  // Cast vectors into arrays such that gsl can cope with them
  const double* sqrts = &sqrts_vec[0];
  const double* photon_momentum = &photon_momentum_vec[0];
  const double* dsigma_dk_pipi_pipi_opp = &dsigma_dk_pipi_pipi_opp_vec[0];
  const double* dsigma_dk_pipi_pipi_same = &dsigma_dk_pipi_pipi_same_vec[0];
  const double* dsigma_dk_pipi0_pipi0 = &dsigma_dk_pipi0_pipi0_vec[0];
  const double* dsigma_dk_pipi_pi0pi0 = &dsigma_dk_pipi_pi0pi0_vec[0];
  const double* dsigma_dk_pi0pi0_pipi = &dsigma_dk_pi0pi0_pipi_vec[0];

  // Get dimensions
  int N = sqrts_vec.size();
  int M = photon_momentum_vec.size();

  // Create accelerator objects (interpolation lookups)
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  // Initialize bicubic spline interpolation
  gsl_spline2d *spline_pipi_pipi_opp = gsl_spline2d_alloc(gsl_interp2d_bicubic, M, N);
  gsl_spline2d *spline_pipi_pipi_same = gsl_spline2d_alloc(gsl_interp2d_bicubic, M, N);
  gsl_spline2d *spline_pipi0_pipi0 = gsl_spline2d_alloc(gsl_interp2d_bicubic, M, N);
  gsl_spline2d *spline_pipi_pi0pi0 = gsl_spline2d_alloc(gsl_interp2d_bicubic, M, N);
  gsl_spline2d *spline_pi0pi0_pipi = gsl_spline2d_alloc(gsl_interp2d_bicubic, M, N);
  gsl_spline2d_init(spline_pipi_pipi_opp, photon_momentum, sqrts, dsigma_dk_pipi_pipi_opp, M, N);
  gsl_spline2d_init(spline_pipi_pipi_same, photon_momentum, sqrts, dsigma_dk_pipi_pipi_same, M, N);
  gsl_spline2d_init(spline_pipi0_pipi0, photon_momentum, sqrts, dsigma_dk_pipi0_pipi0, M, N);
  gsl_spline2d_init(spline_pipi_pi0pi0, photon_momentum, sqrts, dsigma_dk_pipi_pi0pi0, M, N);
  gsl_spline2d_init(spline_pi0pi0_pipi, photon_momentum, sqrts, dsigma_dk_pi0pi0_pipi, M, N);

  double k_val, sqrts_val;
  double dSigmadk_val_pipi_pipi_opp, dSigmadk_val_pipi_pipi_same, dSigmadk_val_pipi0_pipi0, dSigmadk_val_pipi_pi0pi0, dSigmadk_val_pi0pi0_pipi;

  std::ofstream Brems_dSigmadk;
  Brems_dSigmadk.open(dir + "/Brems_dSigmadk_Xsections.txt");

  Brems_dSigmadk << "# k, pi+- pi-+ -> pi+- pi-+ gamma, pi+- pi+- -> pi+- pi+- gamma, pi pi0 -> pi pi0 gamma, pi pi -> pi0 pi0 gamma, pi0 pi0 -> pi pi gamma"
                                << std::endl;
  sqrts_val = 1.0;  // Fix sqrts for one-dimensional plotting
  Brems_dSigmadk << "# sqrts = " << sqrts_val << " GeV"<< std::endl;
  Brems_dSigmadk << "# k in GeV, dSigma/dk in mbarn/GeV" << std::endl;

  for (int i = 0; i < 400; i += 1) {

    k_val = 0.001 + i * 0.001;
    dSigmadk_val_pipi_pipi_opp = gsl_spline2d_eval(spline_pipi_pipi_opp, k_val, sqrts_val, xacc, yacc);
    dSigmadk_val_pipi_pipi_same = gsl_spline2d_eval(spline_pipi_pipi_same, k_val, sqrts_val, xacc, yacc);
    dSigmadk_val_pipi0_pipi0 = gsl_spline2d_eval(spline_pipi0_pipi0, k_val, sqrts_val, xacc, yacc);
    dSigmadk_val_pipi_pi0pi0 = gsl_spline2d_eval(spline_pipi_pi0pi0, k_val, sqrts_val, xacc, yacc);
    dSigmadk_val_pi0pi0_pipi = gsl_spline2d_eval(spline_pi0pi0_pipi, k_val, sqrts_val, xacc, yacc);

    Brems_dSigmadk << k_val << "\t" << dSigmadk_val_pipi_pipi_opp
        << "\t" << dSigmadk_val_pipi_pipi_same << "\t" << dSigmadk_val_pipi0_pipi0
        << "\t" << dSigmadk_val_pipi_pi0pi0 << "\t" << dSigmadk_val_pi0pi0_pipi << std::endl;
  }
  Brems_dSigmadk.close();

  gsl_spline2d_free(spline_pipi_pipi_opp);
  gsl_spline2d_free(spline_pipi_pipi_same);
  gsl_spline2d_free(spline_pipi0_pipi0);
  gsl_spline2d_free(spline_pipi_pi0pi0);
  gsl_spline2d_free(spline_pi0pi0_pipi);
  gsl_interp_accel_free (xacc);
  gsl_interp_accel_free (yacc);
}

void Brems_writing::write_dSigma_dTheta(std::string dir) {
  // Read in tabularized differential cross sections as vectors
  std::vector<double> sqrts_vec = BREMS_SQRTS;
  std::vector<double> photon_angle_vec = BREMS_THETA;
  std::vector<double> dsigma_dtheta_pipi_pipi_opp_vec = BREMS_PIPI_PIPI_OPP_DIFF_SIG_THETA;
  std::vector<double> dsigma_dtheta_pipi_pipi_same_vec = BREMS_PIPI_PIPI_SAME_DIFF_SIG_THETA;
  std::vector<double> dsigma_dtheta_pipi0_pipi0_vec = BREMS_PIPI0_PIPI0_DIFF_SIG_THETA;
  std::vector<double> dsigma_dtheta_pipi_pi0pi0_vec = BREMS_PIPI_PI0PI0_DIFF_SIG_THETA;
  std::vector<double> dsigma_dtheta_pi0pi0_pipi_vec = BREMS_PI0PI0_PIPI_DIFF_SIG_THETA;

  // Cast vectors into arrays such that gsl can cope with them
  const double* sqrts = &sqrts_vec[0];
  const double* photon_angle = &photon_angle_vec[0];
  const double* dsigma_dtheta_pipi_pipi_opp = &dsigma_dtheta_pipi_pipi_opp_vec[0];
  const double* dsigma_dtheta_pipi_pipi_same = &dsigma_dtheta_pipi_pipi_same_vec[0];
  const double* dsigma_dtheta_pipi0_pipi0 = &dsigma_dtheta_pipi0_pipi0_vec[0];
  const double* dsigma_dtheta_pipi_pi0pi0 = &dsigma_dtheta_pipi_pi0pi0_vec[0];
  const double* dsigma_dtheta_pi0pi0_pipi = &dsigma_dtheta_pi0pi0_pipi_vec[0];

  // Get dimensions
  int N = sqrts_vec.size();
  int M = photon_angle_vec.size();

  // Create accelerator objects (interpolation lookups)
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  // Initialize bicubic spline interpolation
  gsl_spline2d *spline_pipi_pipi_opp = gsl_spline2d_alloc(gsl_interp2d_bicubic, M, N);
  gsl_spline2d *spline_pipi_pipi_same = gsl_spline2d_alloc(gsl_interp2d_bicubic, M, N);
  gsl_spline2d *spline_pipi0_pipi0 = gsl_spline2d_alloc(gsl_interp2d_bicubic, M, N);
  gsl_spline2d *spline_pipi_pi0pi0 = gsl_spline2d_alloc(gsl_interp2d_bicubic, M, N);
  gsl_spline2d *spline_pi0pi0_pipi = gsl_spline2d_alloc(gsl_interp2d_bicubic, M, N);
  gsl_spline2d_init(spline_pipi_pipi_opp, photon_angle, sqrts, dsigma_dtheta_pipi_pipi_opp, M, N);
  gsl_spline2d_init(spline_pipi_pipi_same, photon_angle, sqrts, dsigma_dtheta_pipi_pipi_same, M, N);
  gsl_spline2d_init(spline_pipi0_pipi0, photon_angle, sqrts, dsigma_dtheta_pipi0_pipi0, M, N);
  gsl_spline2d_init(spline_pipi_pi0pi0, photon_angle, sqrts, dsigma_dtheta_pipi_pi0pi0, M, N);
  gsl_spline2d_init(spline_pi0pi0_pipi, photon_angle, sqrts, dsigma_dtheta_pi0pi0_pipi, M, N);

  double theta_val, sqrts_val;
  double dSigmadtheta_val_pipi_pipi_opp, dSigmadtheta_val_pipi_pipi_same, dSigmadtheta_val_pipi0_pipi0, dSigmadtheta_val_pipi_pi0pi0, dSigmadtheta_val_pi0pi0_pipi;
  sqrts_val = 1.0;  // Fix sqrts for one-dimensional plotting

  std::ofstream Brems_dSigmadTheta;
  Brems_dSigmadTheta.open(dir + "/Brems_dSigmadTheta_Xsections.txt");
  Brems_dSigmadTheta << "# theta, pi+- pi-+ -> pi+- pi-+ gamma, pi+- pi+- -> pi+- pi+- gamma, pi pi0 -> pi pi0 gamma, pi pi -> pi0 pi0 gamma, pi0 pi0 -> pi pi gamma"
                                << std::endl;
  Brems_dSigmadTheta << "# sqrts = " << sqrts_val << " GeV"<< std::endl;
  Brems_dSigmadTheta << "# dSigma/dk in mbarn" << std::endl;

  double step = M_PI / 100.0;
  for (int i = 0; i < 100; i += 1) {
    theta_val = 0.0 + i * step;
    dSigmadtheta_val_pipi_pipi_opp = gsl_spline2d_eval(spline_pipi_pipi_opp, theta_val, sqrts_val, xacc, yacc);
    dSigmadtheta_val_pipi_pipi_same = gsl_spline2d_eval(spline_pipi_pipi_same, theta_val, sqrts_val, xacc, yacc);
    dSigmadtheta_val_pipi0_pipi0 = gsl_spline2d_eval(spline_pipi0_pipi0, theta_val, sqrts_val, xacc, yacc);
    dSigmadtheta_val_pipi_pi0pi0 = gsl_spline2d_eval(spline_pipi_pi0pi0, theta_val, sqrts_val, xacc, yacc);
    dSigmadtheta_val_pi0pi0_pipi = gsl_spline2d_eval(spline_pi0pi0_pipi, theta_val, sqrts_val, xacc, yacc);

    Brems_dSigmadTheta << theta_val << "\t" << dSigmadtheta_val_pipi_pipi_opp
        << "\t" << dSigmadtheta_val_pipi_pipi_same << "\t" << dSigmadtheta_val_pipi0_pipi0
        << "\t" << dSigmadtheta_val_pipi_pi0pi0 << "\t" << dSigmadtheta_val_pi0pi0_pipi << std::endl;
  }
  Brems_dSigmadTheta.close();

  gsl_spline2d_free(spline_pipi_pipi_opp);
  gsl_spline2d_free(spline_pipi_pipi_same);
  gsl_spline2d_free(spline_pipi0_pipi0);
  gsl_spline2d_free(spline_pipi_pi0pi0);
  gsl_spline2d_free(spline_pi0pi0_pipi);
  gsl_interp_accel_free (xacc);
  gsl_interp_accel_free (yacc);
}
