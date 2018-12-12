// common_functions.h

#ifndef COMMON_FUNCTIONS_H_
#define COMMON_FUNCTIONS_H_

#include <cmath>
#include <fstream>
#include <iostream>

/* Heaviside step function, necessary for processes with different kinematic
 * thresholds for broad rho mesons. Sometimes the s-channel is not accessible
 * in a certain energy region while the t- and u- channel is. Thus, the
 * s-channels need to be excluded in the low-energy region.
 */

inline double HeavisideTheta(double x) {
  if (x >= 0.0) {
    return 1.0;
  } else {
    return 0.0;
  }
}

/* Kinematics:
 * m1: mass of the first incoming particle
 * m2: mass of the second incoming particle
 * m3: mass of the first outgoing particle
 * s: squared center-of-mass energy of the binary scattering process
 * The photon as the second outgoing particle is not considered, as it is
 * massless.
 */

// Minimum value of Mandelstam t
inline double min_mandelstam_t(double s, double m1, double m2, double m3) {
  double tmin;
  tmin = pow(m1, 2) + pow(m3, 2) -
         2 * (s + pow(m1, 2) - pow(m2, 2)) * (s + pow(m3, 2)) / (4 * s) -
         2 * sqrt(pow((s + pow(m1, 2) - pow(m2, 2)), 2) - 4 * s * pow(m1, 2)) *
             (s - pow(m3, 2)) / (4 * s);
  return tmin;
}

// Maximum value of Mandelstam t
inline double max_mandelstam_t(double s, double m1, double m2, double m3) {
  double tmax;
  tmax = pow(m1, 2) + pow(m3, 2) -
         2 * (s + pow(m1, 2) - pow(m2, 2)) * (s + pow(m3, 2)) / (4 * s) +
         2 * sqrt(pow((s + pow(m1, 2) - pow(m2, 2)), 2) - 4 * s * pow(m1, 2)) *
             (s - pow(m3, 2)) / (4 * s);
  return tmax;
}

/* Scattering Angle:
 * Determine the scattering angle theta (in degrees) from kinematic properties.
 * m1: mass of the first incoming particle
 * m2: mass of the second incoming particle
 * m3: mass of the first outgoing particle
 * sqrts: center-of-mass energy
 * pcm: center-of-mass momentum
 */

inline double theta_from_t(double t, double sqrts, double m1, double m2,
                           double m3) {
  double s = sqrts * sqrts;
  double pcm = sqrt(pow(pow(m1, 2) - pow(m2, 2), 2) -
                    2 * (pow(m1, 2) + pow(m2, 2)) * s + pow(s, 2)) /
               (2. * sqrt(s));
  double costheta =
      (t - pow(m2, 2) +
       0.5 * (s + pow(m2, 2) - pow(m1, 2)) * (s - pow(m3, 2)) / s) /
      (pcm * (s - pow(m3, 2)) / sqrts);

  double theta = acos(costheta) * 180.0 / M_PI;

  return theta;
}

// Constant parameters:
const double Const = 0.059;
const double g_POR = 11.93; // 22.6 in case no form factors are applied
const double ma1 = 1.26;
const double ghat = 6.4483;
const double eta1 = 2.22388;
const double eta2 = 2.39014;
const double delta = -0.6426;
const double C4 = -0.14095;
const double Gammaa1 = 0.4;
const double momega = 0.783;
const double mpion = 0.138;
const double Pi = M_PI;
const double mrho = 0.776; // Rho meson assumed to be stable (pole mass)

#endif
