/*
 * Collection of total cross sections as a function of the squared
 * center-of-mass energy (s). The cross section does in principle also depend
 * on the mass of the rho meson, it's value is defined in "common_functions.h".
 * The naming conventions (C11, C12, ...) follow from the PhD thesis of
 * S. Turbide (See README for link).
 */

#include "include/total_Xsections.h"
#include "include/common_functions.h"

//C11
double total_xsection_C11(double s) {
  double sigma;
  double spin_deg_factor = 3.0;
  double tmin = min_mandelstam_t(s, mpion, mrho, mpion);
  double tmax = max_mandelstam_t(s, mpion, mrho, mpion);

  sigma =
      (pow(Const, 2) * pow(ghat, 4) *
       ((pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(ma1, 8) + pow(mpion, 8) - pow(mpion, 4) * pow(mrho, 4) -
               2 * pow(ma1, 2) * pow(mpion, 2) *
                   (pow(mpion, 2) - pow(mrho, 2)) * (pow(mrho, 2) + s) +
               pow(ma1, 6) * (-4 * pow(mpion, 2) + 2 * s) +
               pow(ma1, 4) * (4 * pow(mpion, 4) - pow(mrho, 4) +
                              2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) -
                              2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(ma1, 8) +
               pow(mpion, 4) * pow(pow(mpion, 2) - pow(mrho, 2), 2) +
               2 * pow(ma1, 6) * (-2 * pow(mpion, 2) + pow(mrho, 2) + s) +
               2 * pow(ma1, 2) * pow(mpion, 2) *
                   (-pow(mrho, 4) + pow(mpion, 2) * (2 * pow(mrho, 2) - s) +
                    pow(mrho, 2) * s) +
               pow(ma1, 4) *
                   (4 * pow(mpion, 4) + pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                    2 * pow(s, 2) - 4 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
          pow(eta1, 2) *
              (pow(ma1, 8) + pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) -
               2 * pow(ma1, 6) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
               2 * pow(mpion, 2) * pow(mrho, 4) * s +
               pow(mpion, 4) * (3 * pow(mrho, 4) + 2 * pow(mrho, 2) * s) +
               pow(ma1, 4) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                              pow(mpion, 2) * (8 * pow(mrho, 2) - 4 * s) -
                              4 * pow(mrho, 2) * s + 2 * pow(s, 2)) -
               2 * pow(ma1, 2) *
                   (pow(mrho, 2) * s * (-pow(mrho, 2) + s) +
                    pow(mpion, 4) * (3 * pow(mrho, 2) + s) +
                    pow(mpion, 2) *
                        (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s))))) /
            ((pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(ma1, 2) - tmax)) +
        (8 * pow(-2 + delta, 2) * pow(mpion, 2) *
         (4 * pow(mpion, 2) - pow(mrho, 2))) /
            ((pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(mpion, 2) - tmax)) -
        (8 * pow(-2 + delta, 2) * pow(mpion, 2) * tmax) /
            (pow(mrho, 2) * pow(pow(mpion, 2) - s, 2)) -
        (8 * pow(-2 + delta, 2) * pow(mpion, 2) * tmax) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) *
         (-8 * C4 * pow(mrho, 4) +
          pow(mpion, 2) * (2 + delta - 8 * C4 * pow(mrho, 2)) -
          (2 + 3 * delta) * s + pow(mrho, 2) * (-2 + 3 * delta + 16 * C4 * s)) *
         tmax) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(mpion, 2) + pow(mrho, 2) - 2 * s) * (pow(mpion, 2) + s) +
          eta1 * (-2 * pow(mpion, 4) + pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                  2 * pow(s, 2) + pow(mpion, 2) * (pow(mrho, 2) + s))) *
         tmax) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) *
              (3 * pow(ma1, 4) + 4 * pow(mpion, 4) + pow(mrho, 4) +
               pow(mpion, 2) * (8 * pow(mrho, 2) - 4 * s) -
               4 * pow(ma1, 2) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
               4 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
          pow(eta2, 2) *
              (3 * pow(ma1, 4) + 4 * pow(mpion, 4) + pow(mrho, 4) -
               2 * pow(mrho, 2) * s + 2 * pow(s, 2) -
               4 * pow(mpion, 2) * (pow(mrho, 2) + s) +
               4 * pow(ma1, 2) * (-2 * pow(mpion, 2) + pow(mrho, 2) + s)) -
          2 * eta1 * eta2 *
              (3 * pow(ma1, 4) + 4 * pow(mpion, 4) - pow(mrho, 4) +
               2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) -
               2 * pow(mrho, 2) * s + 2 * pow(s, 2) +
               pow(ma1, 2) * (-8 * pow(mpion, 2) + 4 * s))) *
         tmax) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) +
        (8 *
         (pow(delta, 2) * (8 * pow(mpion, 4) + 3 * pow(mrho, 4) +
                           4 * pow(mpion, 2) * (3 * pow(mrho, 2) - 2 * s) -
                           6 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
          4 * pow(mrho, 4) *
              (3 + 12 * C4 * (2 * pow(mpion, 2) - s) +
               8 * pow(C4, 2) * pow(-2 * pow(mpion, 2) + s, 2)) -
          4 * delta * pow(mrho, 2) *
              (16 * C4 * pow(mpion, 4) +
               2 * pow(mpion, 2) * (3 + 6 * C4 * pow(mrho, 2) - 8 * C4 * s) +
               pow(mrho, 2) * (3 - 6 * C4 * s) + s * (-3 + 4 * C4 * s))) *
         tmax) /
            (pow(mrho, 4) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) *
         (pow(mpion, 4) * (-2 + 3 * delta - 8 * C4 * pow(mrho, 2)) +
          (pow(mrho, 2) - s) * ((-2 + 3 * delta) * s +
                                pow(mrho, 2) * (-2 + delta - 8 * C4 * s)) +
          4 * pow(mpion, 2) *
              (2 * C4 * pow(mrho, 4) + delta * s -
               pow(mrho, 2) * (-1 + delta + 4 * C4 * s))) *
         tmax) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(mpion, 2) + s) *
              (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
               (pow(mrho, 2) - s) * s) +
          eta1 * (-4 * pow(mpion, 6) + pow(pow(mrho, 2) - s, 2) * s +
                  pow(mpion, 4) * (3 * pow(mrho, 2) + s) -
                  pow(mpion, 2) *
                      (pow(mrho, 4) - pow(mrho, 2) * s + 2 * pow(s, 2)))) *
         tmax) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(mpion, 8) - pow(mrho, 4) * pow(s, 2) + pow(s, 4) -
               pow(mpion, 4) *
                   (pow(mrho, 4) + 2 * pow(mrho, 2) * s - 4 * pow(s, 2)) +
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * s - 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(s, 2) * pow(pow(mrho, 2) + s, 2) +
               pow(mpion, 4) * pow(pow(mrho, 2) + 2 * s, 2) -
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + 2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta1, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) -
               4 * pow(mpion, 2) * pow(pow(mrho, 2) - s, 2) * s +
               pow(pow(mrho, 2) - s, 2) * pow(s, 2) +
               pow(mpion, 4) *
                   (3 * pow(mrho, 4) - 6 * pow(mrho, 2) * s + 4 * pow(s, 2)))) *
         tmax) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta1, 2) *
              (pow(ma1, 4) * s + pow(mpion, 4) * (-3 * pow(mrho, 2) + 2 * s) +
               s * (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s + pow(s, 2)) -
               2 * pow(mpion, 2) *
                   (pow(mrho, 4) - 4 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
               pow(ma1, 2) * (2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
                              3 * s * (-pow(mrho, 2) + s))) -
          2 * eta1 * eta2 *
              (pow(ma1, 4) * s +
               s * (2 * pow(mpion, 4) + 4 * pow(mpion, 2) * (pow(mrho, 2) - s) +
                    s * (-2 * pow(mrho, 2) + s)) +
               pow(ma1, 2) * (pow(mpion, 2) * (pow(mrho, 2) - 4 * s) +
                              s * (-2 * pow(mrho, 2) + 3 * s))) +
          pow(eta2, 2) *
              (-4 * pow(mpion, 2) * s * (pow(ma1, 2) + pow(mrho, 2) + s) +
               pow(mpion, 4) * (pow(mrho, 2) + 2 * s) +
               s * (pow(ma1, 4) + s * (pow(mrho, 2) + s) +
                    pow(ma1, 2) * (pow(mrho, 2) + 3 * s)))) *
         tmax) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta1 * (-4 * pow(mpion, 4) *
                      (6 * C4 * pow(mrho, 4) + 2 * delta * s +
                       pow(mrho, 2) * (1 - 2 * delta - 8 * C4 * s)) +
                  2 * pow(mpion, 2) *
                      (4 * delta * pow(s, 2) +
                       pow(mrho, 2) * s * (6 - 7 * delta - 16 * C4 * s) +
                       2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * s)) -
                  (pow(mrho, 2) - s) * s *
                      (-2 * delta * s +
                       pow(mrho, 2) * (-6 + 3 * delta + 8 * C4 * s))) +
          eta2 * (delta * (2 * pow(mpion, 4) * (pow(mrho, 2) + 4 * s) +
                           pow(mpion, 2) * (2 * pow(mrho, 4) +
                                            pow(mrho, 2) * s - 8 * pow(s, 2)) +
                           s * (-2 * pow(mrho, 4) - pow(mrho, 2) * s +
                                2 * pow(s, 2))) -
                  2 * pow(mrho, 2) *
                      (4 * C4 * pow(mpion, 4) * (pow(mrho, 2) + 4 * s) +
                       pow(mpion, 2) * (s * (5 - 16 * C4 * s) +
                                        pow(mrho, 2) * (2 - 8 * C4 * s)) +
                       s * (s * (-3 + 4 * C4 * s) +
                            pow(mrho, 2) * (-2 + 4 * C4 * s))))) *
         tmax) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (eta1 - eta2) *
         (delta *
              (eta1 *
                   (4 * pow(mpion, 6) +
                    pow(mpion, 4) * (7 * pow(mrho, 2) - 8 * s) +
                    pow(ma1, 4) * (pow(mpion, 2) - s) -
                    pow(ma1, 2) * (2 * pow(mpion, 2) + pow(mrho, 2) - 2 * s) *
                        (2 * pow(mpion, 2) - s) +
                    pow(mpion, 2) * s * (-8 * pow(mrho, 2) + 5 * s) +
                    s * (pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2))) +
               eta2 *
                   (-4 * pow(mpion, 6) -
                    pow(mpion, 4) * (pow(mrho, 2) - 8 * s) +
                    pow(ma1, 4) * (-pow(mpion, 2) + s) +
                    pow(mpion, 2) * (2 * pow(mrho, 4) - 5 * pow(s, 2)) +
                    s * (-2 * pow(mrho, 4) + pow(mrho, 2) * s + pow(s, 2)) +
                    pow(ma1, 2) * (4 * pow(mpion, 4) - 6 * pow(mpion, 2) * s +
                                   s * (pow(mrho, 2) + 2 * s)))) -
          2 * pow(mrho, 2) *
              (eta1 *
                   (8 * C4 * pow(mpion, 6) +
                    pow(mpion, 4) * (3 + 8 * C4 * (pow(mrho, 2) - 2 * s)) +
                    2 * C4 * pow(ma1, 4) * (pow(mpion, 2) - s) +
                    2 * pow(mpion, 2) * s *
                        (-1 - 6 * C4 * pow(mrho, 2) + 5 * C4 * s) -
                    pow(ma1, 2) *
                        (8 * C4 * pow(mpion, 4) +
                         pow(mpion, 2) * (1 + 2 * C4 * (pow(mrho, 2) - 6 * s)) +
                         2 * C4 * s * (-pow(mrho, 2) + 2 * s)) +
                    s * (-(s * (1 + 2 * C4 * s)) +
                         pow(mrho, 2) * (1 + 4 * C4 * s))) +
               eta2 * (2 * C4 * pow(ma1, 4) * (-pow(mpion, 2) + s) -
                       (pow(mpion, 2) - s) *
                           (8 * C4 * pow(mpion, 4) - 2 * pow(mrho, 2) + s +
                            2 * C4 * pow(s, 2) +
                            pow(mpion, 2) *
                                (3 - 4 * C4 * (pow(mrho, 2) + 2 * s))) +
                       pow(ma1, 2) *
                           (8 * C4 * pow(mpion, 4) +
                            2 * C4 * s * (pow(mrho, 2) + 2 * s) +
                            pow(mpion, 2) *
                                (1 - 2 * C4 * (pow(mrho, 2) + 6 * s)))))) *
         tmax) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) * (delta - 4 * C4 * pow(mrho, 2)) * pow(tmax, 2)) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (16 * (-2 + delta) * (delta - 4 * C4 * pow(mrho, 2)) * s *
         pow(tmax, 2)) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) * (pow(mrho, 2) - s) + 2 * eta1 * eta2 * s -
          pow(eta2, 2) * s) *
         (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
          (pow(mrho, 2) - s) * s) *
         pow(tmax, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (-(eta1 * (pow(mpion, 2) + 2 * pow(mrho, 2) - 3 * s)) -
          eta2 * (pow(mpion, 2) + s)) *
         pow(tmax, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (pow(mpion, 2) + s) *
         (-2 * eta2 * s + eta1 * (pow(mpion, 2) - pow(mrho, 2) + s)) *
         pow(tmax, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 3) *
         (eta1 * (pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s) -
          eta2 * (pow(ma1, 2) - 2 * pow(mpion, 2) + pow(mrho, 2) + s)) *
         pow(tmax, 2)) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) -
        (8 * (delta - 4 * C4 * pow(mrho, 2)) *
         (delta * (4 * pow(mpion, 2) + 3 * pow(mrho, 2) - 2 * s) -
          2 * pow(mrho, 2) * (3 + 8 * C4 * pow(mpion, 2) - 4 * C4 * s)) *
         pow(tmax, 2)) /
            (pow(mrho, 4) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta2, 2) * s *
              (pow(ma1, 2) - 4 * pow(mpion, 2) + pow(mrho, 2) + 3 * s) +
          pow(eta1, 2) * (2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
                          s * (pow(ma1, 2) - 3 * pow(mrho, 2) + 3 * s)) -
          2 * eta1 * eta2 *
              (pow(mpion, 2) * (pow(mrho, 2) - 4 * s) +
               s * (pow(ma1, 2) - 2 * pow(mrho, 2) + 3 * s))) *
         pow(tmax, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (2 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta1 * (4 * delta * pow(s, 2) -
                  2 * pow(mrho, 2) * s * (-2 + 3 * delta + 8 * C4 * s) +
                  pow(mrho, 4) * (-2 + delta + 16 * C4 * s) -
                  2 * pow(mpion, 2) *
                      (8 * C4 * pow(mrho, 4) + 4 * delta * s +
                       pow(mrho, 2) * (2 - 3 * delta - 16 * C4 * s))) +
          eta2 * (pow(mpion, 2) * (8 * delta * s +
                                   pow(mrho, 2) * (-2 + delta - 32 * C4 * s)) +
                  s * (-4 * delta * s +
                       pow(mrho, 2) * (-2 + delta + 16 * C4 * s)))) *
         pow(tmax, 2)) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) *
         (delta * (eta1 * (pow(ma1, 2) * (pow(mpion, 2) - s) -
                           (2 * pow(mpion, 2) + pow(mrho, 2) - 2 * s) *
                               (2 * pow(mpion, 2) - s)) +
                   eta2 * (4 * pow(mpion, 4) - 6 * pow(mpion, 2) * s +
                           pow(ma1, 2) * (-pow(mpion, 2) + s) +
                           s * (pow(mrho, 2) + 2 * s))) +
          2 * pow(mrho, 2) *
              (eta1 *
                   (8 * C4 * pow(mpion, 4) +
                    2 * C4 * s * (pow(ma1, 2) - pow(mrho, 2) + 2 * s) +
                    pow(mpion, 2) *
                        (1 - 2 * C4 * (pow(ma1, 2) - pow(mrho, 2) + 6 * s))) -
               eta2 * (8 * C4 * pow(mpion, 4) +
                       2 * C4 * s * (pow(ma1, 2) + pow(mrho, 2) + 2 * s) -
                       pow(mpion, 2) *
                           (-1 +
                            2 * C4 * (pow(ma1, 2) + pow(mrho, 2) + 6 * s))))) *
         pow(tmax, 2)) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 4) * pow(tmax, 3)) /
            (3. * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                   2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * pow(eta1 - eta2, 2) * (delta - 4 * C4 * pow(mrho, 2)) *
         pow(tmax, 3)) /
            (3. * pow(mrho, 2) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (16 * pow(delta - 4 * C4 * pow(mrho, 2), 2) * pow(tmax, 3)) /
            (3. * pow(mrho, 4) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (4 * (-2 + delta) * eta1 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         pow(tmax, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (2 * pow(eta1 - eta2, 4) * (pow(ma1, 2) - s) * s * pow(tmax, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (2 * pow(eta1 - eta2, 2) * s *
         (-2 * eta1 * eta2 * s + pow(eta2, 2) * s +
          pow(eta1, 2) * (-pow(mrho, 2) + s)) *
         pow(tmax, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (2 * eta2 * (delta - 4 * C4 * pow(mrho, 2)) * s +
          eta1 * (-2 * delta * s + pow(mrho, 2) * (-2 + delta + 8 * C4 * s))) *
         pow(tmax, 3)) /
            (3. * pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(ma1, 8) + pow(mpion, 8) - pow(mpion, 4) * pow(mrho, 4) -
               2 * pow(ma1, 2) * pow(mpion, 2) *
                   (pow(mpion, 2) - pow(mrho, 2)) * (pow(mrho, 2) + s) +
               pow(ma1, 6) * (-4 * pow(mpion, 2) + 2 * s) +
               pow(ma1, 4) * (4 * pow(mpion, 4) - pow(mrho, 4) +
                              2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) -
                              2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(ma1, 8) +
               pow(mpion, 4) * pow(pow(mpion, 2) - pow(mrho, 2), 2) +
               2 * pow(ma1, 6) * (-2 * pow(mpion, 2) + pow(mrho, 2) + s) +
               2 * pow(ma1, 2) * pow(mpion, 2) *
                   (-pow(mrho, 4) + pow(mpion, 2) * (2 * pow(mrho, 2) - s) +
                    pow(mrho, 2) * s) +
               pow(ma1, 4) *
                   (4 * pow(mpion, 4) + pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                    2 * pow(s, 2) - 4 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
          pow(eta1, 2) *
              (pow(ma1, 8) + pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) -
               2 * pow(ma1, 6) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
               2 * pow(mpion, 2) * pow(mrho, 4) * s +
               pow(mpion, 4) * (3 * pow(mrho, 4) + 2 * pow(mrho, 2) * s) +
               pow(ma1, 4) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                              pow(mpion, 2) * (8 * pow(mrho, 2) - 4 * s) -
                              4 * pow(mrho, 2) * s + 2 * pow(s, 2)) -
               2 * pow(ma1, 2) *
                   (pow(mrho, 2) * s * (-pow(mrho, 2) + s) +
                    pow(mpion, 4) * (3 * pow(mrho, 2) + s) +
                    pow(mpion, 2) *
                        (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s))))) /
            ((pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(ma1, 2) - tmin)) -
        (8 * pow(-2 + delta, 2) * pow(mpion, 2) *
         (4 * pow(mpion, 2) - pow(mrho, 2))) /
            ((pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(mpion, 2) - tmin)) +
        (8 * pow(-2 + delta, 2) * pow(mpion, 2) * tmin) /
            (pow(mrho, 2) * pow(pow(mpion, 2) - s, 2)) +
        (8 * pow(-2 + delta, 2) * pow(mpion, 2) * tmin) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) *
         (-8 * C4 * pow(mrho, 4) +
          pow(mpion, 2) * (2 + delta - 8 * C4 * pow(mrho, 2)) -
          (2 + 3 * delta) * s + pow(mrho, 2) * (-2 + 3 * delta + 16 * C4 * s)) *
         tmin) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(mpion, 2) + pow(mrho, 2) - 2 * s) * (pow(mpion, 2) + s) +
          eta1 * (-2 * pow(mpion, 4) + pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                  2 * pow(s, 2) + pow(mpion, 2) * (pow(mrho, 2) + s))) *
         tmin) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) *
              (3 * pow(ma1, 4) + 4 * pow(mpion, 4) + pow(mrho, 4) +
               pow(mpion, 2) * (8 * pow(mrho, 2) - 4 * s) -
               4 * pow(ma1, 2) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) -
               4 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
          pow(eta2, 2) *
              (3 * pow(ma1, 4) + 4 * pow(mpion, 4) + pow(mrho, 4) -
               2 * pow(mrho, 2) * s + 2 * pow(s, 2) -
               4 * pow(mpion, 2) * (pow(mrho, 2) + s) +
               4 * pow(ma1, 2) * (-2 * pow(mpion, 2) + pow(mrho, 2) + s)) -
          2 * eta1 * eta2 *
              (3 * pow(ma1, 4) + 4 * pow(mpion, 4) - pow(mrho, 4) +
               2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) -
               2 * pow(mrho, 2) * s + 2 * pow(s, 2) +
               pow(ma1, 2) * (-8 * pow(mpion, 2) + 4 * s))) *
         tmin) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) -
        (8 *
         (pow(delta, 2) * (8 * pow(mpion, 4) + 3 * pow(mrho, 4) +
                           4 * pow(mpion, 2) * (3 * pow(mrho, 2) - 2 * s) -
                           6 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
          4 * pow(mrho, 4) *
              (3 + 12 * C4 * (2 * pow(mpion, 2) - s) +
               8 * pow(C4, 2) * pow(-2 * pow(mpion, 2) + s, 2)) -
          4 * delta * pow(mrho, 2) *
              (16 * C4 * pow(mpion, 4) +
               2 * pow(mpion, 2) * (3 + 6 * C4 * pow(mrho, 2) - 8 * C4 * s) +
               pow(mrho, 2) * (3 - 6 * C4 * s) + s * (-3 + 4 * C4 * s))) *
         tmin) /
            (pow(mrho, 4) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) *
         (pow(mpion, 4) * (-2 + 3 * delta - 8 * C4 * pow(mrho, 2)) +
          (pow(mrho, 2) - s) * ((-2 + 3 * delta) * s +
                                pow(mrho, 2) * (-2 + delta - 8 * C4 * s)) +
          4 * pow(mpion, 2) *
              (2 * C4 * pow(mrho, 4) + delta * s -
               pow(mrho, 2) * (-1 + delta + 4 * C4 * s))) *
         tmin) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(mpion, 2) + s) *
              (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
               (pow(mrho, 2) - s) * s) +
          eta1 * (-4 * pow(mpion, 6) + pow(pow(mrho, 2) - s, 2) * s +
                  pow(mpion, 4) * (3 * pow(mrho, 2) + s) -
                  pow(mpion, 2) *
                      (pow(mrho, 4) - pow(mrho, 2) * s + 2 * pow(s, 2)))) *
         tmin) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(mpion, 8) - pow(mrho, 4) * pow(s, 2) + pow(s, 4) -
               pow(mpion, 4) *
                   (pow(mrho, 4) + 2 * pow(mrho, 2) * s - 4 * pow(s, 2)) +
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * s - 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(s, 2) * pow(pow(mrho, 2) + s, 2) +
               pow(mpion, 4) * pow(pow(mrho, 2) + 2 * s, 2) -
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + 2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta1, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) -
               4 * pow(mpion, 2) * pow(pow(mrho, 2) - s, 2) * s +
               pow(pow(mrho, 2) - s, 2) * pow(s, 2) +
               pow(mpion, 4) *
                   (3 * pow(mrho, 4) - 6 * pow(mrho, 2) * s + 4 * pow(s, 2)))) *
         tmin) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta1, 2) *
              (pow(ma1, 4) * s + pow(mpion, 4) * (-3 * pow(mrho, 2) + 2 * s) +
               s * (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s + pow(s, 2)) -
               2 * pow(mpion, 2) *
                   (pow(mrho, 4) - 4 * pow(mrho, 2) * s + 2 * pow(s, 2)) +
               pow(ma1, 2) * (2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
                              3 * s * (-pow(mrho, 2) + s))) -
          2 * eta1 * eta2 *
              (pow(ma1, 4) * s +
               s * (2 * pow(mpion, 4) + 4 * pow(mpion, 2) * (pow(mrho, 2) - s) +
                    s * (-2 * pow(mrho, 2) + s)) +
               pow(ma1, 2) * (pow(mpion, 2) * (pow(mrho, 2) - 4 * s) +
                              s * (-2 * pow(mrho, 2) + 3 * s))) +
          pow(eta2, 2) *
              (-4 * pow(mpion, 2) * s * (pow(ma1, 2) + pow(mrho, 2) + s) +
               pow(mpion, 4) * (pow(mrho, 2) + 2 * s) +
               s * (pow(ma1, 4) + s * (pow(mrho, 2) + s) +
                    pow(ma1, 2) * (pow(mrho, 2) + 3 * s)))) *
         tmin) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta1 * (4 * pow(mpion, 4) *
                      (6 * C4 * pow(mrho, 4) + 2 * delta * s +
                       pow(mrho, 2) * (1 - 2 * delta - 8 * C4 * s)) -
                  2 * pow(mpion, 2) *
                      (4 * delta * pow(s, 2) +
                       pow(mrho, 2) * s * (6 - 7 * delta - 16 * C4 * s) +
                       2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * s)) +
                  (pow(mrho, 2) - s) * s *
                      (-2 * delta * s +
                       pow(mrho, 2) * (-6 + 3 * delta + 8 * C4 * s))) +
          eta2 *
              (-(delta *
                 (2 * pow(mpion, 4) * (pow(mrho, 2) + 4 * s) +
                  pow(mpion, 2) *
                      (2 * pow(mrho, 4) + pow(mrho, 2) * s - 8 * pow(s, 2)) +
                  s * (-2 * pow(mrho, 4) - pow(mrho, 2) * s + 2 * pow(s, 2)))) +
               2 * pow(mrho, 2) *
                   (4 * C4 * pow(mpion, 4) * (pow(mrho, 2) + 4 * s) +
                    pow(mpion, 2) * (s * (5 - 16 * C4 * s) +
                                     pow(mrho, 2) * (2 - 8 * C4 * s)) +
                    s * (s * (-3 + 4 * C4 * s) +
                         pow(mrho, 2) * (-2 + 4 * C4 * s))))) *
         tmin) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * (eta1 - eta2) *
         (delta *
              (eta1 *
                   (4 * pow(mpion, 6) +
                    pow(mpion, 4) * (7 * pow(mrho, 2) - 8 * s) +
                    pow(ma1, 4) * (pow(mpion, 2) - s) -
                    pow(ma1, 2) * (2 * pow(mpion, 2) + pow(mrho, 2) - 2 * s) *
                        (2 * pow(mpion, 2) - s) +
                    pow(mpion, 2) * s * (-8 * pow(mrho, 2) + 5 * s) +
                    s * (pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2))) +
               eta2 *
                   (-4 * pow(mpion, 6) -
                    pow(mpion, 4) * (pow(mrho, 2) - 8 * s) +
                    pow(ma1, 4) * (-pow(mpion, 2) + s) +
                    pow(mpion, 2) * (2 * pow(mrho, 4) - 5 * pow(s, 2)) +
                    s * (-2 * pow(mrho, 4) + pow(mrho, 2) * s + pow(s, 2)) +
                    pow(ma1, 2) * (4 * pow(mpion, 4) - 6 * pow(mpion, 2) * s +
                                   s * (pow(mrho, 2) + 2 * s)))) -
          2 * pow(mrho, 2) *
              (eta1 *
                   (8 * C4 * pow(mpion, 6) +
                    pow(mpion, 4) * (3 + 8 * C4 * (pow(mrho, 2) - 2 * s)) +
                    2 * C4 * pow(ma1, 4) * (pow(mpion, 2) - s) +
                    2 * pow(mpion, 2) * s *
                        (-1 - 6 * C4 * pow(mrho, 2) + 5 * C4 * s) -
                    pow(ma1, 2) *
                        (8 * C4 * pow(mpion, 4) +
                         pow(mpion, 2) * (1 + 2 * C4 * (pow(mrho, 2) - 6 * s)) +
                         2 * C4 * s * (-pow(mrho, 2) + 2 * s)) +
                    s * (-(s * (1 + 2 * C4 * s)) +
                         pow(mrho, 2) * (1 + 4 * C4 * s))) +
               eta2 * (2 * C4 * pow(ma1, 4) * (-pow(mpion, 2) + s) -
                       (pow(mpion, 2) - s) *
                           (8 * C4 * pow(mpion, 4) - 2 * pow(mrho, 2) + s +
                            2 * C4 * pow(s, 2) +
                            pow(mpion, 2) *
                                (3 - 4 * C4 * (pow(mrho, 2) + 2 * s))) +
                       pow(ma1, 2) *
                           (8 * C4 * pow(mpion, 4) +
                            2 * C4 * s * (pow(mrho, 2) + 2 * s) +
                            pow(mpion, 2) *
                                (1 - 2 * C4 * (pow(mrho, 2) + 6 * s)))))) *
         tmin) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) * (delta - 4 * C4 * pow(mrho, 2)) * pow(tmin, 2)) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (16 * (-2 + delta) * (delta - 4 * C4 * pow(mrho, 2)) * s *
         pow(tmin, 2)) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 2) *
         (pow(eta1, 2) * (pow(mrho, 2) - s) + 2 * eta1 * eta2 * s -
          pow(eta2, 2) * s) *
         (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
          (pow(mrho, 2) - s) * s) *
         pow(tmin, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (-(eta1 * (pow(mpion, 2) + 2 * pow(mrho, 2) - 3 * s)) -
          eta2 * (pow(mpion, 2) + s)) *
         pow(tmin, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (2 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (pow(mpion, 2) + s) *
         (-2 * eta2 * s + eta1 * (pow(mpion, 2) - pow(mrho, 2) + s)) *
         pow(tmin, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 3) *
         (-(eta1 * (pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s)) +
          eta2 * (pow(ma1, 2) - 2 * pow(mpion, 2) + pow(mrho, 2) + s)) *
         pow(tmin, 2)) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) +
        (8 * (delta - 4 * C4 * pow(mrho, 2)) *
         (delta * (4 * pow(mpion, 2) + 3 * pow(mrho, 2) - 2 * s) -
          2 * pow(mrho, 2) * (3 + 8 * C4 * pow(mpion, 2) - 4 * C4 * s)) *
         pow(tmin, 2)) /
            (pow(mrho, 4) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta2, 2) * s *
              (pow(ma1, 2) - 4 * pow(mpion, 2) + pow(mrho, 2) + 3 * s) +
          pow(eta1, 2) * (2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
                          s * (pow(ma1, 2) - 3 * pow(mrho, 2) + 3 * s)) -
          2 * eta1 * eta2 *
              (pow(mpion, 2) * (pow(mrho, 2) - 4 * s) +
               s * (pow(ma1, 2) - 2 * pow(mrho, 2) + 3 * s))) *
         pow(tmin, 2)) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (2 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta1 * (4 * delta * pow(s, 2) -
                  2 * pow(mrho, 2) * s * (-2 + 3 * delta + 8 * C4 * s) +
                  pow(mrho, 4) * (-2 + delta + 16 * C4 * s) -
                  2 * pow(mpion, 2) *
                      (8 * C4 * pow(mrho, 4) + 4 * delta * s +
                       pow(mrho, 2) * (2 - 3 * delta - 16 * C4 * s))) +
          eta2 * (pow(mpion, 2) * (8 * delta * s +
                                   pow(mrho, 2) * (-2 + delta - 32 * C4 * s)) +
                  s * (-4 * delta * s +
                       pow(mrho, 2) * (-2 + delta + 16 * C4 * s)))) *
         pow(tmin, 2)) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (4 * (eta1 - eta2) *
         (delta * (eta1 * (pow(ma1, 2) * (pow(mpion, 2) - s) -
                           (2 * pow(mpion, 2) + pow(mrho, 2) - 2 * s) *
                               (2 * pow(mpion, 2) - s)) +
                   eta2 * (4 * pow(mpion, 4) - 6 * pow(mpion, 2) * s +
                           pow(ma1, 2) * (-pow(mpion, 2) + s) +
                           s * (pow(mrho, 2) + 2 * s))) +
          2 * pow(mrho, 2) *
              (eta1 *
                   (8 * C4 * pow(mpion, 4) +
                    2 * C4 * s * (pow(ma1, 2) - pow(mrho, 2) + 2 * s) +
                    pow(mpion, 2) *
                        (1 - 2 * C4 * (pow(ma1, 2) - pow(mrho, 2) + 6 * s))) -
               eta2 * (8 * C4 * pow(mpion, 4) +
                       2 * C4 * s * (pow(ma1, 2) + pow(mrho, 2) + 2 * s) -
                       pow(mpion, 2) *
                           (-1 +
                            2 * C4 * (pow(ma1, 2) + pow(mrho, 2) + 6 * s))))) *
         pow(tmin, 2)) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (pow(eta1 - eta2, 4) * pow(tmin, 3)) /
            (3. * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                   2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * pow(eta1 - eta2, 2) * (delta - 4 * C4 * pow(mrho, 2)) *
         pow(tmin, 3)) /
            (3. * pow(mrho, 2) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (16 * pow(delta - 4 * C4 * pow(mrho, 2), 2) * pow(tmin, 3)) /
            (3. * pow(mrho, 4) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (-2 + delta) * eta1 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         pow(tmin, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (2 * pow(eta1 - eta2, 4) * (pow(ma1, 2) - s) * s * pow(tmin, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (2 * pow(eta1 - eta2, 2) * s *
         (-2 * eta1 * eta2 * s + pow(eta2, 2) * s +
          pow(eta1, 2) * (-pow(mrho, 2) + s)) *
         pow(tmin, 3)) /
            (3. * (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (4 * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (2 * eta2 * (delta - 4 * C4 * pow(mrho, 2)) * s +
          eta1 * (-2 * delta * s + pow(mrho, 2) * (-2 + delta + 8 * C4 * s))) *
         pow(tmin, 3)) /
            (3. * pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (2 * pow(eta1 - eta2, 2) *
         (pow(eta1, 2) *
              (2 * pow(ma1, 6) -
               3 * pow(ma1, 4) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) +
               pow(mrho, 2) * (pow(mrho, 2) - s) * s -
               pow(mpion, 4) * (3 * pow(mrho, 2) + s) +
               pow(mpion, 2) * (-2 * pow(mrho, 4) + 3 * pow(mrho, 2) * s) +
               pow(ma1, 2) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                              pow(mpion, 2) * (8 * pow(mrho, 2) - 4 * s) -
                              4 * pow(mrho, 2) * s + 2 * pow(s, 2))) -
          2 * eta1 * eta2 *
              (2 * pow(ma1, 6) -
               pow(mpion, 2) * (pow(mpion, 2) - pow(mrho, 2)) *
                   (pow(mrho, 2) + s) +
               pow(ma1, 4) * (-6 * pow(mpion, 2) + 3 * s) +
               pow(ma1, 2) * (4 * pow(mpion, 4) - pow(mrho, 4) +
                              2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) -
                              2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (2 * pow(ma1, 6) +
               3 * pow(ma1, 4) * (-2 * pow(mpion, 2) + pow(mrho, 2) + s) +
               pow(mpion, 2) *
                   (-pow(mrho, 4) + pow(mpion, 2) * (2 * pow(mrho, 2) - s) +
                    pow(mrho, 2) * s) +
               pow(ma1, 2) *
                   (4 * pow(mpion, 4) + pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                    2 * pow(s, 2) - 4 * pow(mpion, 2) * (pow(mrho, 2) + s)))) *
         log(abs(-pow(ma1, 2) + tmax))) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) -
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (-2 * eta1 * eta2 *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               2 * pow(ma1, 2) * pow(mpion, 4) * s +
               pow(mpion, 2) * (pow(ma1, 4) * (pow(mrho, 2) - 4 * s) +
                                4 * pow(ma1, 2) * (pow(mrho, 2) - s) * s +
                                pow(mrho, 2) * pow(s, 2)) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + s * (-2 * pow(mrho, 2) + s) +
                    pow(ma1, 2) * (-2 * pow(mrho, 2) + 3 * s))) +
          pow(eta2, 2) *
              (pow(mpion, 8) -
               4 * pow(ma1, 2) * pow(mpion, 2) * s *
                   (pow(ma1, 2) + pow(mrho, 2) + s) +
               pow(mpion, 4) *
                   (pow(mrho, 2) * s + pow(ma1, 2) * (pow(mrho, 2) + 2 * s)) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + s * (pow(mrho, 2) + s) +
                    pow(ma1, 2) * (pow(mrho, 2) + 3 * s))) +
          pow(eta1, 2) *
              (pow(mpion, 8) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + 2 * pow(mrho, 4) -
                    3 * pow(ma1, 2) * (pow(mrho, 2) - s) -
                    3 * pow(mrho, 2) * s + pow(s, 2)) +
               pow(mpion, 4) * (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                                pow(ma1, 2) * (-3 * pow(mrho, 2) + 2 * s)) +
               2 * pow(mpion, 2) *
                   (pow(ma1, 4) * (pow(mrho, 2) - 2 * s) +
                    pow(mrho, 2) * s * (-pow(mrho, 2) + s) -
                    pow(ma1, 2) * (pow(mrho, 4) - 4 * pow(mrho, 2) * s +
                                   2 * pow(s, 2))))) *
         log(abs(-pow(ma1, 2) + tmax))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (eta1 - eta2) *
         (delta *
              (eta2 *
                   (pow(mpion, 6) * pow(mrho, 2) * (2 * pow(mpion, 2) - s) +
                    pow(ma1, 8) * (-pow(mpion, 2) + s) +
                    pow(ma1, 6) * (5 * pow(mpion, 4) - 7 * pow(mpion, 2) * s +
                                   s * (pow(mrho, 2) + 2 * s)) +
                    pow(ma1, 4) *
                        (-8 * pow(mpion, 6) -
                         pow(mpion, 4) * (pow(mrho, 2) - 14 * s) +
                         pow(mpion, 2) * (2 * pow(mrho, 4) - pow(mrho, 2) * s -
                                          7 * pow(s, 2)) +
                         s * (-2 * pow(mrho, 4) + pow(mrho, 2) * s +
                              pow(s, 2))) +
                    pow(ma1, 2) * pow(mpion, 2) *
                        (4 * pow(mpion, 6) +
                         pow(mpion, 4) * (pow(mrho, 2) - 8 * s) +
                         s * (2 * pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2)) +
                         pow(mpion, 2) *
                             (-2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                              5 * pow(s, 2)))) +
               eta1 *
                   (pow(ma1, 8) * (pow(mpion, 2) - s) +
                    pow(ma1, 6) *
                        (-5 * pow(mpion, 4) + (pow(mrho, 2) - 2 * s) * s +
                         pow(mpion, 2) * (-2 * pow(mrho, 2) + 7 * s)) +
                    pow(mpion, 2) * pow(mrho, 2) *
                        (2 * pow(mpion, 6) +
                         pow(mpion, 4) * (4 * pow(mrho, 2) - 5 * s) +
                         pow(mrho, 4) * s -
                         pow(mpion, 2) * (pow(mrho, 4) + 3 * pow(mrho, 2) * s -
                                          2 * pow(s, 2))) +
                    pow(ma1, 4) *
                        (8 * pow(mpion, 6) +
                         pow(mpion, 4) * (9 * pow(mrho, 2) - 14 * s) +
                         pow(mpion, 2) * s * (-9 * pow(mrho, 2) + 7 * s) +
                         s * (pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2))) +
                    pow(ma1, 2) *
                        (-4 * pow(mpion, 8) +
                         pow(mrho, 4) * s * (-pow(mrho, 2) + s) +
                         pow(mpion, 6) * (-11 * pow(mrho, 2) + 8 * s) +
                         pow(mpion, 4) *
                             (-3 * pow(mrho, 4) + 17 * pow(mrho, 2) * s -
                              5 * pow(s, 2)) +
                         pow(mpion, 2) *
                             (pow(mrho, 6) - 5 * pow(mrho, 2) * pow(s, 2) +
                              pow(s, 3))))) -
          2 * pow(mrho, 2) *
              (eta2 *
                   (pow(mpion, 8) * (1 + 2 * C4 * pow(mrho, 2)) -
                    2 * C4 * pow(mpion, 6) * pow(mrho, 2) * s +
                    2 * C4 * pow(ma1, 8) * (-pow(mpion, 2) + s) +
                    pow(ma1, 4) *
                        (-16 * C4 * pow(mpion, 6) +
                         pow(mpion, 4) *
                             (-4 + 6 * C4 * pow(mrho, 2) + 28 * C4 * s) +
                         2 * pow(mpion, 2) *
                             (pow(mrho, 2) + s - 3 * C4 * pow(mrho, 2) * s -
                              7 * C4 * pow(s, 2)) +
                         s * (-2 * pow(mrho, 2) + s + 2 * C4 * pow(s, 2))) +
                    pow(ma1, 6) * (10 * C4 * pow(mpion, 4) +
                                   2 * C4 * s * (pow(mrho, 2) + 2 * s) +
                                   pow(mpion, 2) *
                                       (1 - 2 * C4 * (pow(mrho, 2) + 7 * s))) +
                    pow(ma1, 2) * pow(mpion, 2) *
                        (8 * C4 * pow(mpion, 6) -
                         2 * pow(mpion, 4) *
                             (-2 + 3 * C4 * pow(mrho, 2) + 8 * C4 * s) +
                         s * (2 * pow(mrho, 2) + s - 2 * C4 * pow(s, 2)) +
                         2 * pow(mpion, 2) *
                             (pow(mrho, 2) * (-1 + 3 * C4 * s) +
                              s * (-3 + 5 * C4 * s)))) +
               eta1 *
                   (pow(mpion, 8) * (-1 + 6 * C4 * pow(mrho, 2)) +
                    2 * C4 * pow(ma1, 8) * (pow(mpion, 2) - s) +
                    pow(mpion, 2) * pow(mrho, 4) * s +
                    2 * pow(mpion, 6) * pow(mrho, 2) * (2 - 5 * C4 * s) -
                    pow(ma1, 6) *
                        (10 * C4 * pow(mpion, 4) +
                         pow(mpion, 2) * (1 + 2 * C4 * (pow(mrho, 2) - 7 * s)) +
                         2 * C4 * s * (-pow(mrho, 2) + 2 * s)) -
                    pow(mpion, 4) * pow(mrho, 2) *
                        (pow(mrho, 2) + s * (3 - 4 * C4 * s)) +
                    pow(ma1, 4) *
                        (16 * C4 * pow(mpion, 6) +
                         2 * pow(mpion, 4) *
                             (2 + 5 * C4 * pow(mrho, 2) - 14 * C4 * s) +
                         2 * pow(mpion, 2) * s *
                             (-1 - 7 * C4 * pow(mrho, 2) + 7 * C4 * s) +
                         s * (-(s * (1 + 2 * C4 * s)) +
                              pow(mrho, 2) * (1 + 4 * C4 * s))) -
                    pow(ma1, 2) *
                        (8 * C4 * pow(mpion, 8) +
                         pow(mrho, 2) * (pow(mrho, 2) - s) * s +
                         2 * pow(mpion, 6) *
                             (2 + 7 * C4 * pow(mrho, 2) - 8 * C4 * s) +
                         pow(mpion, 2) * (-pow(mrho, 4) + pow(s, 2) +
                                          8 * C4 * pow(mrho, 2) * pow(s, 2) -
                                          2 * C4 * pow(s, 3)) +
                         pow(mpion, 4) * (pow(mrho, 2) * (3 - 22 * C4 * s) +
                                          2 * s * (-3 + 5 * C4 * s)))))) *
         log(abs(-pow(ma1, 2) + tmax))) /
            ((pow(ma1, 2) - pow(mpion, 2)) * pow(mrho, 2) *
             (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (16 * pow(-2 + delta, 2) * pow(mpion, 2) *
         log(abs(-pow(mpion, 2) + tmax))) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) -
        (8 * pow(-2 + delta, 2) *
         (3 * pow(mpion, 4) - 4 * pow(mpion, 2) * (pow(mrho, 2) - s) +
          pow(pow(mrho, 2) - s, 2)) *
         log(abs(-pow(mpion, 2) + tmax))) /
            ((pow(mpion, 2) - s) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                                    2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) * (eta1 - eta2) * pow(mpion, 2) *
         (2 * eta1 * pow(mpion, 2) - 2 * eta2 * pow(mpion, 2) -
          eta1 * pow(mrho, 2)) *
         (pow(mpion, 2) - s) * log(abs(-pow(mpion, 2) + tmax))) /
            ((pow(ma1, 2) - pow(mpion, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) * (eta1 - eta2) * pow(mpion, 2) * (pow(ma1, 2) - s) *
         (pow(mpion, 2) - s) *
         (-(eta2 * (pow(mpion, 2) + s)) +
          eta1 * (pow(mpion, 2) - pow(mrho, 2) + s)) *
         log(abs(-pow(mpion, 2) + tmax))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) *
         (-(delta * (4 * pow(mpion, 2) - pow(mrho, 2)) *
            (pow(mpion, 2) + pow(mrho, 2) - s)) +
          2 * pow(mrho, 2) *
              (8 * C4 * pow(mpion, 4) - pow(mrho, 2) + s +
               pow(mpion, 2) * (3 - 8 * C4 * s))) *
         log(abs(-pow(mpion, 2) + tmax))) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (2 * pow(eta1 - eta2, 2) *
         (pow(eta1, 2) *
              (2 * pow(ma1, 6) -
               3 * pow(ma1, 4) * (2 * pow(mpion, 2) + pow(mrho, 2) - s) +
               pow(mrho, 2) * (pow(mrho, 2) - s) * s -
               pow(mpion, 4) * (3 * pow(mrho, 2) + s) +
               pow(mpion, 2) * (-2 * pow(mrho, 4) + 3 * pow(mrho, 2) * s) +
               pow(ma1, 2) * (4 * pow(mpion, 4) + pow(mrho, 4) +
                              pow(mpion, 2) * (8 * pow(mrho, 2) - 4 * s) -
                              4 * pow(mrho, 2) * s + 2 * pow(s, 2))) -
          2 * eta1 * eta2 *
              (2 * pow(ma1, 6) -
               pow(mpion, 2) * (pow(mpion, 2) - pow(mrho, 2)) *
                   (pow(mrho, 2) + s) +
               pow(ma1, 4) * (-6 * pow(mpion, 2) + 3 * s) +
               pow(ma1, 2) * (4 * pow(mpion, 4) - pow(mrho, 4) +
                              2 * pow(mpion, 2) * (pow(mrho, 2) - 2 * s) -
                              2 * pow(mrho, 2) * s + 2 * pow(s, 2))) +
          pow(eta2, 2) *
              (2 * pow(ma1, 6) +
               3 * pow(ma1, 4) * (-2 * pow(mpion, 2) + pow(mrho, 2) + s) +
               pow(mpion, 2) *
                   (-pow(mrho, 4) + pow(mpion, 2) * (2 * pow(mrho, 2) - s) +
                    pow(mrho, 2) * s) +
               pow(ma1, 2) *
                   (4 * pow(mpion, 4) + pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                    2 * pow(s, 2) - 4 * pow(mpion, 2) * (pow(mrho, 2) + s)))) *
         log(abs(-pow(ma1, 2) + tmin))) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) +
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (-2 * eta1 * eta2 *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               2 * pow(ma1, 2) * pow(mpion, 4) * s +
               pow(mpion, 2) * (pow(ma1, 4) * (pow(mrho, 2) - 4 * s) +
                                4 * pow(ma1, 2) * (pow(mrho, 2) - s) * s +
                                pow(mrho, 2) * pow(s, 2)) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + s * (-2 * pow(mrho, 2) + s) +
                    pow(ma1, 2) * (-2 * pow(mrho, 2) + 3 * s))) +
          pow(eta2, 2) *
              (pow(mpion, 8) -
               4 * pow(ma1, 2) * pow(mpion, 2) * s *
                   (pow(ma1, 2) + pow(mrho, 2) + s) +
               pow(mpion, 4) *
                   (pow(mrho, 2) * s + pow(ma1, 2) * (pow(mrho, 2) + 2 * s)) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + s * (pow(mrho, 2) + s) +
                    pow(ma1, 2) * (pow(mrho, 2) + 3 * s))) +
          pow(eta1, 2) *
              (pow(mpion, 8) +
               pow(ma1, 2) * s *
                   (pow(ma1, 4) + 2 * pow(mrho, 4) -
                    3 * pow(ma1, 2) * (pow(mrho, 2) - s) -
                    3 * pow(mrho, 2) * s + pow(s, 2)) +
               pow(mpion, 4) * (2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                                pow(ma1, 2) * (-3 * pow(mrho, 2) + 2 * s)) +
               2 * pow(mpion, 2) *
                   (pow(ma1, 4) * (pow(mrho, 2) - 2 * s) +
                    pow(mrho, 2) * s * (-pow(mrho, 2) + s) -
                    pow(ma1, 2) * (pow(mrho, 4) - 4 * pow(mrho, 2) * s +
                                   2 * pow(s, 2))))) *
         log(abs(-pow(ma1, 2) + tmin))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * (eta1 - eta2) *
         (delta *
              (eta2 *
                   (pow(mpion, 6) * pow(mrho, 2) * (2 * pow(mpion, 2) - s) +
                    pow(ma1, 8) * (-pow(mpion, 2) + s) +
                    pow(ma1, 6) * (5 * pow(mpion, 4) - 7 * pow(mpion, 2) * s +
                                   s * (pow(mrho, 2) + 2 * s)) +
                    pow(ma1, 4) *
                        (-8 * pow(mpion, 6) -
                         pow(mpion, 4) * (pow(mrho, 2) - 14 * s) +
                         pow(mpion, 2) * (2 * pow(mrho, 4) - pow(mrho, 2) * s -
                                          7 * pow(s, 2)) +
                         s * (-2 * pow(mrho, 4) + pow(mrho, 2) * s +
                              pow(s, 2))) +
                    pow(ma1, 2) * pow(mpion, 2) *
                        (4 * pow(mpion, 6) +
                         pow(mpion, 4) * (pow(mrho, 2) - 8 * s) +
                         s * (2 * pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2)) +
                         pow(mpion, 2) *
                             (-2 * pow(mrho, 4) - 3 * pow(mrho, 2) * s +
                              5 * pow(s, 2)))) +
               eta1 *
                   (pow(ma1, 8) * (pow(mpion, 2) - s) +
                    pow(ma1, 6) *
                        (-5 * pow(mpion, 4) + (pow(mrho, 2) - 2 * s) * s +
                         pow(mpion, 2) * (-2 * pow(mrho, 2) + 7 * s)) +
                    pow(mpion, 2) * pow(mrho, 2) *
                        (2 * pow(mpion, 6) +
                         pow(mpion, 4) * (4 * pow(mrho, 2) - 5 * s) +
                         pow(mrho, 4) * s -
                         pow(mpion, 2) * (pow(mrho, 4) + 3 * pow(mrho, 2) * s -
                                          2 * pow(s, 2))) +
                    pow(ma1, 4) *
                        (8 * pow(mpion, 6) +
                         pow(mpion, 4) * (9 * pow(mrho, 2) - 14 * s) +
                         pow(mpion, 2) * s * (-9 * pow(mrho, 2) + 7 * s) +
                         s * (pow(mrho, 4) + pow(mrho, 2) * s - pow(s, 2))) +
                    pow(ma1, 2) *
                        (-4 * pow(mpion, 8) +
                         pow(mrho, 4) * s * (-pow(mrho, 2) + s) +
                         pow(mpion, 6) * (-11 * pow(mrho, 2) + 8 * s) +
                         pow(mpion, 4) *
                             (-3 * pow(mrho, 4) + 17 * pow(mrho, 2) * s -
                              5 * pow(s, 2)) +
                         pow(mpion, 2) *
                             (pow(mrho, 6) - 5 * pow(mrho, 2) * pow(s, 2) +
                              pow(s, 3))))) -
          2 * pow(mrho, 2) *
              (eta2 *
                   (pow(mpion, 8) * (1 + 2 * C4 * pow(mrho, 2)) -
                    2 * C4 * pow(mpion, 6) * pow(mrho, 2) * s +
                    2 * C4 * pow(ma1, 8) * (-pow(mpion, 2) + s) +
                    pow(ma1, 4) *
                        (-16 * C4 * pow(mpion, 6) +
                         pow(mpion, 4) *
                             (-4 + 6 * C4 * pow(mrho, 2) + 28 * C4 * s) +
                         2 * pow(mpion, 2) *
                             (pow(mrho, 2) + s - 3 * C4 * pow(mrho, 2) * s -
                              7 * C4 * pow(s, 2)) +
                         s * (-2 * pow(mrho, 2) + s + 2 * C4 * pow(s, 2))) +
                    pow(ma1, 6) * (10 * C4 * pow(mpion, 4) +
                                   2 * C4 * s * (pow(mrho, 2) + 2 * s) +
                                   pow(mpion, 2) *
                                       (1 - 2 * C4 * (pow(mrho, 2) + 7 * s))) +
                    pow(ma1, 2) * pow(mpion, 2) *
                        (8 * C4 * pow(mpion, 6) -
                         2 * pow(mpion, 4) *
                             (-2 + 3 * C4 * pow(mrho, 2) + 8 * C4 * s) +
                         s * (2 * pow(mrho, 2) + s - 2 * C4 * pow(s, 2)) +
                         2 * pow(mpion, 2) *
                             (pow(mrho, 2) * (-1 + 3 * C4 * s) +
                              s * (-3 + 5 * C4 * s)))) +
               eta1 *
                   (pow(mpion, 8) * (-1 + 6 * C4 * pow(mrho, 2)) +
                    2 * C4 * pow(ma1, 8) * (pow(mpion, 2) - s) +
                    pow(mpion, 2) * pow(mrho, 4) * s +
                    2 * pow(mpion, 6) * pow(mrho, 2) * (2 - 5 * C4 * s) -
                    pow(ma1, 6) *
                        (10 * C4 * pow(mpion, 4) +
                         pow(mpion, 2) * (1 + 2 * C4 * (pow(mrho, 2) - 7 * s)) +
                         2 * C4 * s * (-pow(mrho, 2) + 2 * s)) -
                    pow(mpion, 4) * pow(mrho, 2) *
                        (pow(mrho, 2) + s * (3 - 4 * C4 * s)) +
                    pow(ma1, 4) *
                        (16 * C4 * pow(mpion, 6) +
                         2 * pow(mpion, 4) *
                             (2 + 5 * C4 * pow(mrho, 2) - 14 * C4 * s) +
                         2 * pow(mpion, 2) * s *
                             (-1 - 7 * C4 * pow(mrho, 2) + 7 * C4 * s) +
                         s * (-(s * (1 + 2 * C4 * s)) +
                              pow(mrho, 2) * (1 + 4 * C4 * s))) -
                    pow(ma1, 2) *
                        (8 * C4 * pow(mpion, 8) +
                         pow(mrho, 2) * (pow(mrho, 2) - s) * s +
                         2 * pow(mpion, 6) *
                             (2 + 7 * C4 * pow(mrho, 2) - 8 * C4 * s) +
                         pow(mpion, 2) * (-pow(mrho, 4) + pow(s, 2) +
                                          8 * C4 * pow(mrho, 2) * pow(s, 2) -
                                          2 * C4 * pow(s, 3)) +
                         pow(mpion, 4) * (pow(mrho, 2) * (3 - 22 * C4 * s) +
                                          2 * s * (-3 + 5 * C4 * s)))))) *
         log(abs(-pow(ma1, 2) + tmin))) /
            ((pow(ma1, 2) - pow(mpion, 2)) * pow(mrho, 2) *
             (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (16 * pow(-2 + delta, 2) * pow(mpion, 2) *
         log(abs(-pow(mpion, 2) + tmin))) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)) +
        (8 * pow(-2 + delta, 2) *
         (3 * pow(mpion, 4) - 4 * pow(mpion, 2) * (pow(mrho, 2) - s) +
          pow(pow(mrho, 2) - s, 2)) *
         log(abs(-pow(mpion, 2) + tmin))) /
            ((pow(mpion, 2) - s) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                                    2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) * (eta1 - eta2) * pow(mpion, 2) *
         (2 * eta1 * pow(mpion, 2) - 2 * eta2 * pow(mpion, 2) -
          eta1 * pow(mrho, 2)) *
         (pow(mpion, 2) - s) * log(abs(-pow(mpion, 2) + tmin))) /
            ((pow(ma1, 2) - pow(mpion, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) -
        (8 * (-2 + delta) * (eta1 - eta2) * pow(mpion, 2) * (pow(ma1, 2) - s) *
         (pow(mpion, 2) - s) *
         (-(eta2 * (pow(mpion, 2) + s)) +
          eta1 * (pow(mpion, 2) - pow(mrho, 2) + s)) *
         log(abs(-pow(mpion, 2) + tmin))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (8 * (-2 + delta) *
         (delta * (4 * pow(mpion, 2) - pow(mrho, 2)) *
              (pow(mpion, 2) + pow(mrho, 2) - s) -
          2 * pow(mrho, 2) *
              (8 * C4 * pow(mpion, 4) - pow(mrho, 2) + s +
               pow(mpion, 2) * (3 - 8 * C4 * s))) *
         log(abs(-pow(mpion, 2) + tmin))) /
            (pow(mrho, 2) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))))) /
      (512. * Pi);

  return sigma / spin_deg_factor;
}

//C12
double total_xsection_C12(double s) {
  double sigma;
  double spin_deg_factor = 3.0;
  double tmin = min_mandelstam_t(s, mpion, mrho, mpion);
  double tmax = max_mandelstam_t(s, mpion, mrho, mpion);

  sigma = (pow(Const,2)*pow(ghat,4)*(0. - (0.25*pow(-2 + delta,2)*pow(mpion,2)*
           (pow(mpion,4) + pow(pow(mrho,2) - s,2) - 2*pow(mpion,2)*(pow(mrho,2) + s))*tmax)/(pow(mrho,2)*pow(pow(mpion,2) - s,2)) -
        (0.0625*(0. + 8.*pow(mpion,2)*pow(mrho,4) - 12.*pow(mrho,6) + 4.*pow(mrho,4)*s +
             delta*pow(mrho,2)*(-16.*pow(mpion,4) - 16.*pow(mpion,2)*pow(mrho,2) - 4.*pow(mrho,4) + 16.*pow(mrho,2)*s + 4.*pow(s,2)) +
             pow(delta,2)*(8.*pow(mpion,6) + 9.*pow(mrho,6) + pow(mpion,4)*(4.*pow(mrho,2) - 4.*s) - 13.*pow(mrho,4)*s -
                5.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(-2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 2.*pow(s,2))))*tmax)/
         pow(mrho,6) - (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
           (eta2*(pow(mpion,6) + pow(mpion,2)*pow(s,2) + (pow(mrho,2) - 1.*s)*pow(s,2) + pow(mpion,4)*(-1.*pow(mrho,2) + 3.*s)) +
             eta1*(-4.*pow(mpion,6) + pow(mpion,4)*(3.*pow(mrho,2) + s) +
                pow(mpion,2)*(-1.*pow(mrho,4) + pow(mrho,2)*s - 2.*pow(s,2)) + s*(pow(mrho,4) - 2.*pow(mrho,2)*s + pow(s,2))))*tmax)/
         ((-1.*pow(mpion,2) + s)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) +
        (0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta1,2)*(1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) +
                pow(mpion,2)*s*(-4.*pow(mrho,4) + 8.*pow(mrho,2)*s - 4.*pow(s,2)) +
                pow(s,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s + 1.*pow(s,2)) +
                pow(mpion,4)*(3.*pow(mrho,4) - 6.*pow(mrho,2)*s + 4.*pow(s,2))) +
             pow(eta2,2)*(1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) +
                pow(mpion,2)*s*(-2.*pow(mrho,4) - 4.*pow(mrho,2)*s - 4.*pow(s,2)) +
                pow(s,2)*(1.*pow(mrho,4) + 2.*pow(mrho,2)*s + 1.*pow(s,2)) +
                pow(mpion,4)*(1.*pow(mrho,4) + 4.*pow(mrho,2)*s + 4.*pow(s,2))) +
             eta1*eta2*(-2.*pow(mpion,8) + 2.*pow(mrho,4)*pow(s,2) - 2.*pow(s,4) +
                pow(mpion,4)*(2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 8.*pow(s,2)) +
                pow(mpion,2)*s*(-4.*pow(mrho,4) - 4.*pow(mrho,2)*s + 8.*pow(s,2))))*tmax)/
         (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) +
        (0.5*(0. - 4.*C4*pow(mrho,8) - 0.5*pow(mrho,4)*s + pow(mrho,6)*(2. + 2.*C4*s) +
             pow(delta,2)*(2.*pow(mpion,6) + 2.*pow(mrho,6) - 1.5*pow(mpion,4)*s - 2.375*pow(mrho,4)*s - 0.75*pow(mrho,2)*pow(s,2) +
                0.125*pow(s,3) + pow(mpion,2)*(-1.5*pow(mrho,4) + 0.5*pow(mrho,2)*s)) +
             delta*pow(mrho,2)*(-2.*pow(mpion,4) + pow(mpion,2)*(1.*s + pow(mrho,2)*(-1. - 2.*C4*s)) +
                pow(mrho,2)*(2.*C4*pow(mrho,4) + pow(mrho,2)*(-3. + 1.*C4*s) + s*(2. + 1.*C4*s))))*tmax)/pow(mrho,6) -
        (0.5*(pow(mrho,6)*(-1.5 + C4*(-12.*pow(mpion,2) + 6.*s) + pow(C4,2)*(-16.*pow(mpion,4) + 16.*pow(mpion,2)*s - 4.*pow(s,2))) +
             pow(delta,2)*(1.*pow(mpion,6) + 0.125*pow(mrho,6) + pow(mpion,4)*(-2.*pow(mrho,2) - 1.*s) + 0.25*pow(mrho,4)*s -
                0.625*pow(mrho,2)*pow(s,2) + pow(mpion,2)*(-2.5*pow(mrho,4) + 1.75*pow(mrho,2)*s + 0.25*pow(s,2))) +
             delta*pow(mrho,2)*(pow(mpion,4)*(1. + 8.*C4*pow(mrho,2)) +
                pow(mpion,2)*(6.*C4*pow(mrho,4) - 0.5*s + pow(mrho,2)*(3. - 10.*C4*s)) +
                pow(mrho,2)*(pow(mrho,2)*(1.5 - 1.*C4*s) + s*(-2.5 + 3.*C4*s))))*tmax)/pow(mrho,6) -
        (0.25*(pow(delta,2)*(1.*pow(mpion,6) - 1.*pow(mrho,6) + pow(mpion,4)*(-2.499999999999999*pow(mrho,2) - 2.5*s) -
                1.5*pow(mrho,4)*s + 2.*pow(mrho,2)*pow(s,2) - 0.5*pow(s,3) +
                pow(mpion,2)*(3.5*pow(mrho,4) - 1.5000000000000004*pow(mrho,2)*s + 2.*pow(s,2))) +
             pow(mrho,2)*(pow(mpion,4)*(-6. - 8.*C4*pow(mrho,2)) + 2.*pow(s,2) + pow(mrho,4)*(-4. - 8.*C4*s) +
                pow(mrho,2)*s*(-2. + 8.*C4*s) + pow(mpion,2)*(8.*C4*pow(mrho,4) + 4.*s + pow(mrho,2)*(10. - 16.*C4*s))) +
             delta*(-2.*pow(mpion,6) - 5.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) +
                pow(mpion,4)*(8.*pow(mrho,2) + 4.*C4*pow(mrho,4) + 5.*s) + pow(mrho,4)*s*(4. - 4.*C4*s) + pow(mrho,6)*(4. + 4.*C4*s) +
                pow(mpion,2)*(-4.*C4*pow(mrho,6) + 1.*pow(mrho,2)*s - 4.*pow(s,2) + pow(mrho,4)*(-12. + 8.*C4*s))))*tmax)/
         (pow(mrho,4)*(pow(mpion,2) - 1.*s)) + (0.0625*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
           (pow(mrho,2)*(eta2*(4.*pow(mpion,4) - 6.*pow(mpion,2)*s + s*(8.*pow(mrho,2) + 6.*s)) +
                eta1*(-12.*pow(mpion,4) + 4.*pow(mrho,4) + 2.*pow(mrho,2)*s - 6.*pow(s,2) + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s))) +
             delta*(eta1*(8.*pow(mpion,6) - 2.*pow(mrho,6) + pow(mpion,4)*(2.*pow(mrho,2) - 2.*s) - 3.*pow(mrho,4)*s +
                   4.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(2.*pow(mrho,4) - 2.*pow(s,2))) +
                eta2*(pow(mpion,4)*(-2.*pow(mrho,2) - 4.*s) + pow(mpion,2)*s*(3.*pow(mrho,2) + 3.*s) +
                   s*(-4.*pow(mrho,4) - 7.*pow(mrho,2)*s - 1.*pow(s,2)))))*tmax)/
         (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
        (0.1875*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(delta*
              (eta1*(2.6666666666666665*pow(mpion,6) + pow(mpion,4)*(-4.*pow(mrho,2) + 2.*s) +
                   pow(mpion,2)*(-1.3333333333333333*pow(mrho,4) + 6.*pow(mrho,2)*s - 3.3333333333333335*pow(s,2)) +
                   s*(0.3333333333333333*pow(mrho,4) - 1.3333333333333333*pow(mrho,2)*s + 1.*pow(s,2))) +
                eta2*(pow(mpion,4)*(-0.6666666666666666*pow(mrho,2) - 4.*s) +
                   s*(0.6666666666666666*pow(mrho,4) - 1.*pow(mrho,2)*s - 1.*pow(s,2)) +
                   pow(mpion,2)*(-0.6666666666666666*pow(mrho,4) - 0.3333333333333333*pow(mrho,2)*s + 3.6666666666666665*pow(s,2)))) +
             pow(mrho,2)*(eta2*(C4*pow(mpion,4)*(2.6666666666666665*pow(mrho,2) + 10.666666666666666*s) +
                   pow(mpion,2)*(s*(3.3333333333333335 - 10.666666666666666*C4*s) + pow(mrho,2)*(1.3333333333333333 - 5.333333333333333*C4*s)) +
                   s*(s*(-2. + 2.6666666666666665*C4*s) + pow(mrho,2)*(-1.3333333333333333 + 2.6666666666666665*C4*s))) +
                eta1*(pow(mpion,4)*(1.3333333333333333 + 8.*C4*pow(mrho,2) - 10.666666666666666*C4*s) +
                   s*(s*(2. - 2.6666666666666665*C4*s) + pow(mrho,2)*(-2. + 2.6666666666666665*C4*s)) +
                   pow(mpion,2)*(pow(mrho,2)*(2.6666666666666665 - 10.666666666666666*C4*s) + s*(-4. + 10.666666666666666*C4*s)))))*tmax)/
         (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
        (0.0625*(-2. + delta)*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(pow(mpion,2) + s)*
           (-2.*eta2*s + eta1*(pow(mpion,2) - 1.*pow(mrho,2) + s))*pow(tmax,2))/
         ((-1.*pow(mpion,2) + s)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) +
        (0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta1,2)*(pow(mrho,2) - 1.*s) + 2.*eta1*eta2*s - 1.*pow(eta2,2)*s)*
           (pow(mpion,4) + (pow(mrho,2) - 1.*s)*s + pow(mpion,2)*(-1.*pow(mrho,2) + 2.*s))*pow(tmax,2))/
         (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) -
        (0.125*(-1.*pow(mrho,4) + 4.*C4*pow(mrho,6) + delta*pow(mrho,2)*
              (2.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) - 2.*s) +
             pow(delta,2)*(1.*pow(mpion,4) + 0.25*pow(mrho,4) - 1.25*pow(s,2) + pow(mpion,2)*(-3.*pow(mrho,2) + 2.*s)))*pow(tmax,2))/
         pow(mrho,6) + (0.03125*(0. - 4.*pow(mrho,4) + delta*(16.*pow(mrho,4) - 8.*pow(mrho,2)*s) +
             pow(delta,2)*(4.*pow(mpion,4) - 3.*pow(mrho,4) + 2.*pow(mrho,2)*s - 3.*pow(s,2) + pow(mpion,2)*(-4.*pow(mrho,2) + 4.*s)))*
           pow(tmax,2))/pow(mrho,6) + (0.0625*(-32.*C4*pow(mrho,4)*s +
             pow(delta,2)*(1.*pow(mpion,4) + pow(mpion,2)*(-1.0000000000000009*pow(mrho,2) - 2.*s) + s*(-3.*pow(mrho,2) + 1.*s)) +
             delta*(-2.*pow(mpion,4) + (6.*pow(mrho,2) + 16.*C4*pow(mrho,4) - 2.*s)*s + pow(mpion,2)*(2.*pow(mrho,2) + 4.*s)))*
           pow(tmax,2))/(pow(mrho,4)*(pow(mpion,2) - 1.*s)) -
        (0.5625*(C4*pow(mrho,6)*(2.6666666666666665 + 7.111111111111112*C4*pow(mpion,2) - 3.555555555555556*C4*s) +
             pow(delta,2)*(0.11111111111111112*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 0.22222222222222224*s) -
                0.22222222222222224*pow(mrho,2)*s + 0.11111111111111112*pow(s,2)) +
             delta*pow(mrho,2)*(-2.2222222222222223*C4*pow(mrho,4) +
                pow(mpion,2)*(-0.6666666666666666 - 2.6666666666666665*C4*pow(mrho,2)) + 0.22222222222222224*s +
                pow(mrho,2)*(-0.22222222222222224 + 1.777777777777778*C4*s)))*pow(tmax,2))/pow(mrho,6) +
        (0.03125*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(pow(mrho,2)*
              (-2.*eta2*pow(mpion,2) - 5.999999999999999*eta1*pow(mrho,2) + 8.*eta1*s - 2.*eta2*s) +
             delta*(eta1*(-5.999999999999999*pow(mpion,4) + 5.*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) -
                   5.999999999999999*pow(mrho,2)*s + 1.*pow(s,2)) +
                eta2*(4.*pow(mpion,4) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*s) + s*(5.*pow(mrho,2) + 2.*s))))*pow(tmax,2))/
         (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
        (0.15625*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(delta*
              (eta1*(-1.2*pow(mpion,4) + 0.6*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 2.4*s) - 1.6*pow(mrho,2)*s + 1.*pow(s,2)) +
                eta2*(0.8*pow(mpion,4) + (1.*pow(mrho,2) - 0.4*s)*s + pow(mpion,2)*(0.2*pow(mrho,2) + 1.2*s))) +
             pow(mrho,2)*(eta2*(pow(mpion,2)*(-0.4 - 6.4*C4*s) + s*(-0.4 + 3.2*C4*s)) +
                eta1*(s*(0.8 - 3.2*C4*s) + pow(mrho,2)*(-0.4 + 3.2*C4*s) + pow(mpion,2)*(-0.8 - 3.2*C4*pow(mrho,2) + 6.4*C4*s))))*pow(tmax,2)
           )/(pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
        (0.20833333333333331*delta*(-0.8*pow(mrho,2) + 0.8*C4*pow(mrho,4) + delta*(0.8*pow(mpion,2) + 1.*pow(mrho,2) - 0.7*s))*
           pow(tmax,3))/pow(mrho,6) + (0.125*(5.333333333333333*pow(C4,2)*pow(mrho,6) +
             delta*(-0.6666666666666666*pow(mrho,2) - 1.3333333333333333*C4*pow(mrho,4)) +
             pow(delta,2)*(1.*pow(mpion,2) + 1.1666666666666667*pow(mrho,2) - 0.6666666666666666*s))*pow(tmax,3))/pow(mrho,6) +
        (0.10416666666666666*delta*(-0.8*pow(mrho,2) + delta*(0.4*pow(mpion,2) + 1.*pow(mrho,2) - 0.6*s))*pow(tmax,3))/pow(mrho,6) +
        (0.020833333333333332*pow(eta1 - 1.*eta2,2)*s*(-2.*eta1*eta2*s + pow(eta2,2)*s + pow(eta1,2)*(-1.*pow(mrho,2) + s))*pow(tmax,3))/
         (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) +
        (0.10416666666666666*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
           (0.4*eta1*pow(mrho,2) + delta*(-0.2*eta2*pow(mpion,2) - 0.2*eta2*s + eta1*(-0.4*pow(mpion,2) - 0.8*pow(mrho,2) + 1.*s)))*
           pow(tmax,3))/(pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
        (0.14583333333333331*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
           (delta*(-0.14285714285714285*eta2*pow(mpion,2) - 0.42857142857142855*eta2*s +
                eta1*(-0.2857142857142857*pow(mpion,2) - 0.5714285714285714*pow(mrho,2) + 1.*s)) +
             pow(mrho,2)*(1.1428571428571428*C4*eta2*s + eta1*(0.2857142857142857 - 1.1428571428571428*C4*s)))*pow(tmax,3))/
         (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) +
        (0. - 0.5000000000000001*pow(2. - 1.*delta,2)*pow(mpion,4)*pow(mrho,6) + 0.25*pow(2. - 1.*delta,2)*pow(mrho,10) -
           0.5000000000000001*pow(2. - 1.*delta,2)*pow(mrho,6)*pow(s,2) +
           pow(2. - 1.*delta,2)*pow(mpion,2)*pow(mrho,6)*(-1.*pow(mrho,2) + 1.*s))/(pow(mrho,6)*(-2.*pow(mpion,2) + 1.*s + 1.*tmax)) +
        (2.*(0. - 2.*pow(mpion,4)*pow(mrho,4) - 0.5*pow(mrho,8) +
             delta*pow(mrho,4)*(2.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(-2.*pow(mrho,2) - 1.9999999999999998*s)) +
             pow(mpion,2)*(2.*pow(mrho,6) + 2.*pow(mrho,4)*s) +
             pow(delta,2)*pow(mrho,2)*(-2.220446049250313e-16*pow(mpion,6) - 0.125*pow(mrho,6) +
                pow(mpion,4)*(-0.5*pow(mrho,2) + 2.220446049250313e-16*s) + pow(mpion,2)*(0.5*pow(mrho,4) + 0.5*pow(mrho,2)*s)))*
           log(abs(-1.*pow(mpion,2) + 0.5*s + 0.5*tmax)))/(pow(mrho,4)*(pow(mpion,2) - 1.*s)) -
        (0.25*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(eta2*((-2. + 1.*delta)*pow(mpion,6) + (6. - 3.*delta)*pow(mpion,4)*s +
                pow(s,2)*((4. - 2.*delta)*pow(mrho,2) + (2. - 1.*delta)*s) +
                pow(mpion,2)*s*((-4. + 2.*delta)*pow(mrho,2) + (-6. + 3.*delta)*s)) +
             eta1*((2. - 1.*delta)*pow(mpion,6) + (2. - 1.*delta)*pow(mrho,4)*s + (-2. + 1.*delta)*pow(s,3) +
                pow(mpion,4)*((4. - 2.*delta)*pow(mrho,2) + (-6. + 3.*delta)*s) +
                pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (-4. + 2.*delta)*pow(mrho,2)*s + (6. - 3.*delta)*pow(s,2))))*
           log(abs(-2.*pow(mpion,2) + s + tmax)))/(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) +
        (0.25*(0. + 8.*pow(mpion,2)*pow(mrho,6) + 4.*pow(mrho,8) - 8.*pow(mrho,6)*s +
             delta*pow(mrho,4)*(8.*pow(mpion,4) - 8.*pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 16.*s) + 8.*pow(mrho,2)*s +
                8.*pow(s,2)) + pow(delta,2)*pow(mrho,4)*
              (-4.*pow(mpion,4) + 3.*pow(mrho,4) - 2.*pow(mrho,2)*s - 4.*pow(s,2) + pow(mpion,2)*(-6.*pow(mrho,2) + 8.*s)))*
           log(abs(-2.*pow(mpion,2) + 1.*s + 1.*tmax)))/pow(mrho,6) +
        (0.5*(0. + pow(mpion,2)*(4.*pow(mrho,6) - 8.*C4*pow(mrho,8)) - 4.*pow(mrho,6)*s + 8.*C4*pow(mrho,8)*s +
             pow(delta,2)*pow(mrho,4)*(2.*pow(mpion,4) - 1.*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) + 2.*pow(s,2)) +
             delta*pow(mrho,4)*(-4.*pow(mpion,4) + 2.*pow(mrho,2)*s - 4.*pow(s,2) +
                pow(mpion,2)*(-10.*pow(mrho,2) + 4.*C4*pow(mrho,4) + 8.*s) + pow(mrho,4)*(2. - 4.*C4*s)))*
           log(abs(-2.*pow(mpion,2) + 1.*s + 1.*tmax)))/pow(mrho,6)))/
    (16.*Pi*(pow(mpion,4) + pow(pow(mrho,2) - s,2) - 2*pow(mpion,2)*(pow(mrho,2) + s))) -
   (pow(Const,2)*pow(ghat,4)*(0. - (0.25*pow(-2 + delta,2)*pow(mpion,2)*
           (pow(mpion,4) + pow(pow(mrho,2) - s,2) - 2*pow(mpion,2)*(pow(mrho,2) + s))*tmin)/(pow(mrho,2)*pow(pow(mpion,2) - s,2)) -
        (0.0625*(0. + 8.*pow(mpion,2)*pow(mrho,4) - 12.*pow(mrho,6) + 4.*pow(mrho,4)*s +
             delta*pow(mrho,2)*(-16.*pow(mpion,4) - 16.*pow(mpion,2)*pow(mrho,2) - 4.*pow(mrho,4) + 16.*pow(mrho,2)*s + 4.*pow(s,2)) +
             pow(delta,2)*(8.*pow(mpion,6) + 9.*pow(mrho,6) + pow(mpion,4)*(4.*pow(mrho,2) - 4.*s) - 13.*pow(mrho,4)*s -
                5.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(-2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 2.*pow(s,2))))*tmin)/
         pow(mrho,6) - (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
           (eta2*(pow(mpion,6) + pow(mpion,2)*pow(s,2) + (pow(mrho,2) - 1.*s)*pow(s,2) + pow(mpion,4)*(-1.*pow(mrho,2) + 3.*s)) +
             eta1*(-4.*pow(mpion,6) + pow(mpion,4)*(3.*pow(mrho,2) + s) +
                pow(mpion,2)*(-1.*pow(mrho,4) + pow(mrho,2)*s - 2.*pow(s,2)) + s*(pow(mrho,4) - 2.*pow(mrho,2)*s + pow(s,2))))*tmin)/
         ((-1.*pow(mpion,2) + s)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) +
        (0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta1,2)*(1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) +
                pow(mpion,2)*s*(-4.*pow(mrho,4) + 8.*pow(mrho,2)*s - 4.*pow(s,2)) +
                pow(s,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s + 1.*pow(s,2)) +
                pow(mpion,4)*(3.*pow(mrho,4) - 6.*pow(mrho,2)*s + 4.*pow(s,2))) +
             pow(eta2,2)*(1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) +
                pow(mpion,2)*s*(-2.*pow(mrho,4) - 4.*pow(mrho,2)*s - 4.*pow(s,2)) +
                pow(s,2)*(1.*pow(mrho,4) + 2.*pow(mrho,2)*s + 1.*pow(s,2)) +
                pow(mpion,4)*(1.*pow(mrho,4) + 4.*pow(mrho,2)*s + 4.*pow(s,2))) +
             eta1*eta2*(-2.*pow(mpion,8) + 2.*pow(mrho,4)*pow(s,2) - 2.*pow(s,4) +
                pow(mpion,4)*(2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 8.*pow(s,2)) +
                pow(mpion,2)*s*(-4.*pow(mrho,4) - 4.*pow(mrho,2)*s + 8.*pow(s,2))))*tmin)/
         (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) +
        (0.5*(0. - 4.*C4*pow(mrho,8) - 0.5*pow(mrho,4)*s + pow(mrho,6)*(2. + 2.*C4*s) +
             pow(delta,2)*(2.*pow(mpion,6) + 2.*pow(mrho,6) - 1.5*pow(mpion,4)*s - 2.375*pow(mrho,4)*s - 0.75*pow(mrho,2)*pow(s,2) +
                0.125*pow(s,3) + pow(mpion,2)*(-1.5*pow(mrho,4) + 0.5*pow(mrho,2)*s)) +
             delta*pow(mrho,2)*(-2.*pow(mpion,4) + pow(mpion,2)*(1.*s + pow(mrho,2)*(-1. - 2.*C4*s)) +
                pow(mrho,2)*(2.*C4*pow(mrho,4) + pow(mrho,2)*(-3. + 1.*C4*s) + s*(2. + 1.*C4*s))))*tmin)/pow(mrho,6) -
        (0.5*(pow(mrho,6)*(-1.5 + C4*(-12.*pow(mpion,2) + 6.*s) + pow(C4,2)*(-16.*pow(mpion,4) + 16.*pow(mpion,2)*s - 4.*pow(s,2))) +
             pow(delta,2)*(1.*pow(mpion,6) + 0.125*pow(mrho,6) + pow(mpion,4)*(-2.*pow(mrho,2) - 1.*s) + 0.25*pow(mrho,4)*s -
                0.625*pow(mrho,2)*pow(s,2) + pow(mpion,2)*(-2.5*pow(mrho,4) + 1.75*pow(mrho,2)*s + 0.25*pow(s,2))) +
             delta*pow(mrho,2)*(pow(mpion,4)*(1. + 8.*C4*pow(mrho,2)) +
                pow(mpion,2)*(6.*C4*pow(mrho,4) - 0.5*s + pow(mrho,2)*(3. - 10.*C4*s)) +
                pow(mrho,2)*(pow(mrho,2)*(1.5 - 1.*C4*s) + s*(-2.5 + 3.*C4*s))))*tmin)/pow(mrho,6) -
        (0.25*(pow(delta,2)*(1.*pow(mpion,6) - 1.*pow(mrho,6) + pow(mpion,4)*(-2.499999999999999*pow(mrho,2) - 2.5*s) -
                1.5*pow(mrho,4)*s + 2.*pow(mrho,2)*pow(s,2) - 0.5*pow(s,3) +
                pow(mpion,2)*(3.5*pow(mrho,4) - 1.5000000000000004*pow(mrho,2)*s + 2.*pow(s,2))) +
             pow(mrho,2)*(pow(mpion,4)*(-6. - 8.*C4*pow(mrho,2)) + 2.*pow(s,2) + pow(mrho,4)*(-4. - 8.*C4*s) +
                pow(mrho,2)*s*(-2. + 8.*C4*s) + pow(mpion,2)*(8.*C4*pow(mrho,4) + 4.*s + pow(mrho,2)*(10. - 16.*C4*s))) +
             delta*(-2.*pow(mpion,6) - 5.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) +
                pow(mpion,4)*(8.*pow(mrho,2) + 4.*C4*pow(mrho,4) + 5.*s) + pow(mrho,4)*s*(4. - 4.*C4*s) + pow(mrho,6)*(4. + 4.*C4*s) +
                pow(mpion,2)*(-4.*C4*pow(mrho,6) + 1.*pow(mrho,2)*s - 4.*pow(s,2) + pow(mrho,4)*(-12. + 8.*C4*s))))*tmin)/
         (pow(mrho,4)*(pow(mpion,2) - 1.*s)) + (0.0625*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
           (pow(mrho,2)*(eta2*(4.*pow(mpion,4) - 6.*pow(mpion,2)*s + s*(8.*pow(mrho,2) + 6.*s)) +
                eta1*(-12.*pow(mpion,4) + 4.*pow(mrho,4) + 2.*pow(mrho,2)*s - 6.*pow(s,2) + pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s))) +
             delta*(eta1*(8.*pow(mpion,6) - 2.*pow(mrho,6) + pow(mpion,4)*(2.*pow(mrho,2) - 2.*s) - 3.*pow(mrho,4)*s +
                   4.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) + pow(mpion,2)*(2.*pow(mrho,4) - 2.*pow(s,2))) +
                eta2*(pow(mpion,4)*(-2.*pow(mrho,2) - 4.*s) + pow(mpion,2)*s*(3.*pow(mrho,2) + 3.*s) +
                   s*(-4.*pow(mrho,4) - 7.*pow(mrho,2)*s - 1.*pow(s,2)))))*tmin)/
         (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
        (0.1875*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(delta*
              (eta1*(2.6666666666666665*pow(mpion,6) + pow(mpion,4)*(-4.*pow(mrho,2) + 2.*s) +
                   pow(mpion,2)*(-1.3333333333333333*pow(mrho,4) + 6.*pow(mrho,2)*s - 3.3333333333333335*pow(s,2)) +
                   s*(0.3333333333333333*pow(mrho,4) - 1.3333333333333333*pow(mrho,2)*s + 1.*pow(s,2))) +
                eta2*(pow(mpion,4)*(-0.6666666666666666*pow(mrho,2) - 4.*s) +
                   s*(0.6666666666666666*pow(mrho,4) - 1.*pow(mrho,2)*s - 1.*pow(s,2)) +
                   pow(mpion,2)*(-0.6666666666666666*pow(mrho,4) - 0.3333333333333333*pow(mrho,2)*s + 3.6666666666666665*pow(s,2)))) +
             pow(mrho,2)*(eta2*(C4*pow(mpion,4)*(2.6666666666666665*pow(mrho,2) + 10.666666666666666*s) +
                   pow(mpion,2)*(s*(3.3333333333333335 - 10.666666666666666*C4*s) + pow(mrho,2)*(1.3333333333333333 - 5.333333333333333*C4*s)) +
                   s*(s*(-2. + 2.6666666666666665*C4*s) + pow(mrho,2)*(-1.3333333333333333 + 2.6666666666666665*C4*s))) +
                eta1*(pow(mpion,4)*(1.3333333333333333 + 8.*C4*pow(mrho,2) - 10.666666666666666*C4*s) +
                   s*(s*(2. - 2.6666666666666665*C4*s) + pow(mrho,2)*(-2. + 2.6666666666666665*C4*s)) +
                   pow(mpion,2)*(pow(mrho,2)*(2.6666666666666665 - 10.666666666666666*C4*s) + s*(-4. + 10.666666666666666*C4*s)))))*tmin)/
         (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
        (0.0625*(-2. + delta)*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(pow(mpion,2) + s)*
           (-2.*eta2*s + eta1*(pow(mpion,2) - 1.*pow(mrho,2) + s))*pow(tmin,2))/
         ((-1.*pow(mpion,2) + s)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) +
        (0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta1,2)*(pow(mrho,2) - 1.*s) + 2.*eta1*eta2*s - 1.*pow(eta2,2)*s)*
           (pow(mpion,4) + (pow(mrho,2) - 1.*s)*s + pow(mpion,2)*(-1.*pow(mrho,2) + 2.*s))*pow(tmin,2))/
         (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) -
        (0.125*(-1.*pow(mrho,4) + 4.*C4*pow(mrho,6) + delta*pow(mrho,2)*
              (2.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) - 2.*s) +
             pow(delta,2)*(1.*pow(mpion,4) + 0.25*pow(mrho,4) - 1.25*pow(s,2) + pow(mpion,2)*(-3.*pow(mrho,2) + 2.*s)))*pow(tmin,2))/
         pow(mrho,6) + (0.03125*(0. - 4.*pow(mrho,4) + delta*(16.*pow(mrho,4) - 8.*pow(mrho,2)*s) +
             pow(delta,2)*(4.*pow(mpion,4) - 3.*pow(mrho,4) + 2.*pow(mrho,2)*s - 3.*pow(s,2) + pow(mpion,2)*(-4.*pow(mrho,2) + 4.*s)))*
           pow(tmin,2))/pow(mrho,6) + (0.0625*(-32.*C4*pow(mrho,4)*s +
             pow(delta,2)*(1.*pow(mpion,4) + pow(mpion,2)*(-1.0000000000000009*pow(mrho,2) - 2.*s) + s*(-3.*pow(mrho,2) + 1.*s)) +
             delta*(-2.*pow(mpion,4) + (6.*pow(mrho,2) + 16.*C4*pow(mrho,4) - 2.*s)*s + pow(mpion,2)*(2.*pow(mrho,2) + 4.*s)))*
           pow(tmin,2))/(pow(mrho,4)*(pow(mpion,2) - 1.*s)) -
        (0.5625*(C4*pow(mrho,6)*(2.6666666666666665 + 7.111111111111112*C4*pow(mpion,2) - 3.555555555555556*C4*s) +
             pow(delta,2)*(0.11111111111111112*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 0.22222222222222224*s) -
                0.22222222222222224*pow(mrho,2)*s + 0.11111111111111112*pow(s,2)) +
             delta*pow(mrho,2)*(-2.2222222222222223*C4*pow(mrho,4) +
                pow(mpion,2)*(-0.6666666666666666 - 2.6666666666666665*C4*pow(mrho,2)) + 0.22222222222222224*s +
                pow(mrho,2)*(-0.22222222222222224 + 1.777777777777778*C4*s)))*pow(tmin,2))/pow(mrho,6) +
        (0.03125*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(pow(mrho,2)*
              (-2.*eta2*pow(mpion,2) - 5.999999999999999*eta1*pow(mrho,2) + 8.*eta1*s - 2.*eta2*s) +
             delta*(eta1*(-5.999999999999999*pow(mpion,4) + 5.*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) -
                   5.999999999999999*pow(mrho,2)*s + 1.*pow(s,2)) +
                eta2*(4.*pow(mpion,4) + pow(mpion,2)*(1.*pow(mrho,2) - 2.*s) + s*(5.*pow(mrho,2) + 2.*s))))*pow(tmin,2))/
         (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
        (0.15625*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(delta*
              (eta1*(-1.2*pow(mpion,4) + 0.6*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 2.4*s) - 1.6*pow(mrho,2)*s + 1.*pow(s,2)) +
                eta2*(0.8*pow(mpion,4) + (1.*pow(mrho,2) - 0.4*s)*s + pow(mpion,2)*(0.2*pow(mrho,2) + 1.2*s))) +
             pow(mrho,2)*(eta2*(pow(mpion,2)*(-0.4 - 6.4*C4*s) + s*(-0.4 + 3.2*C4*s)) +
                eta1*(s*(0.8 - 3.2*C4*s) + pow(mrho,2)*(-0.4 + 3.2*C4*s) + pow(mpion,2)*(-0.8 - 3.2*C4*pow(mrho,2) + 6.4*C4*s))))*pow(tmin,2)
           )/(pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
        (0.20833333333333331*delta*(-0.8*pow(mrho,2) + 0.8*C4*pow(mrho,4) + delta*(0.8*pow(mpion,2) + 1.*pow(mrho,2) - 0.7*s))*
           pow(tmin,3))/pow(mrho,6) + (0.125*(5.333333333333333*pow(C4,2)*pow(mrho,6) +
             delta*(-0.6666666666666666*pow(mrho,2) - 1.3333333333333333*C4*pow(mrho,4)) +
             pow(delta,2)*(1.*pow(mpion,2) + 1.1666666666666667*pow(mrho,2) - 0.6666666666666666*s))*pow(tmin,3))/pow(mrho,6) +
        (0.10416666666666666*delta*(-0.8*pow(mrho,2) + delta*(0.4*pow(mpion,2) + 1.*pow(mrho,2) - 0.6*s))*pow(tmin,3))/pow(mrho,6) +
        (0.020833333333333332*pow(eta1 - 1.*eta2,2)*s*(-2.*eta1*eta2*s + pow(eta2,2)*s + pow(eta1,2)*(-1.*pow(mrho,2) + s))*pow(tmin,3))/
         (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) +
        (0.10416666666666666*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
           (0.4*eta1*pow(mrho,2) + delta*(-0.2*eta2*pow(mpion,2) - 0.2*eta2*s + eta1*(-0.4*pow(mpion,2) - 0.8*pow(mrho,2) + 1.*s)))*
           pow(tmin,3))/(pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) -
        (0.14583333333333331*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*
           (delta*(-0.14285714285714285*eta2*pow(mpion,2) - 0.42857142857142855*eta2*s +
                eta1*(-0.2857142857142857*pow(mpion,2) - 0.5714285714285714*pow(mrho,2) + 1.*s)) +
             pow(mrho,2)*(1.1428571428571428*C4*eta2*s + eta1*(0.2857142857142857 - 1.1428571428571428*C4*s)))*pow(tmin,3))/
         (pow(mrho,2)*(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2))) +
        (0. - 0.5000000000000001*pow(2. - 1.*delta,2)*pow(mpion,4)*pow(mrho,6) + 0.25*pow(2. - 1.*delta,2)*pow(mrho,10) -
           0.5000000000000001*pow(2. - 1.*delta,2)*pow(mrho,6)*pow(s,2) +
           pow(2. - 1.*delta,2)*pow(mpion,2)*pow(mrho,6)*(-1.*pow(mrho,2) + 1.*s))/(pow(mrho,6)*(-2.*pow(mpion,2) + 1.*s + 1.*tmin)) +
        (2.*(0. - 2.*pow(mpion,4)*pow(mrho,4) - 0.5*pow(mrho,8) +
             delta*pow(mrho,4)*(2.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(-2.*pow(mrho,2) - 1.9999999999999998*s)) +
             pow(mpion,2)*(2.*pow(mrho,6) + 2.*pow(mrho,4)*s) +
             pow(delta,2)*pow(mrho,2)*(-2.220446049250313e-16*pow(mpion,6) - 0.125*pow(mrho,6) +
                pow(mpion,4)*(-0.5*pow(mrho,2) + 2.220446049250313e-16*s) + pow(mpion,2)*(0.5*pow(mrho,4) + 0.5*pow(mrho,2)*s)))*
           log(abs(-1.*pow(mpion,2) + 0.5*s + 0.5*tmin)))/(pow(mrho,4)*(pow(mpion,2) - 1.*s)) -
        (0.25*(eta1 - 1.*eta2)*(pow(ma1,2) - 1.*s)*(eta2*((-2. + 1.*delta)*pow(mpion,6) + (6. - 3.*delta)*pow(mpion,4)*s +
                pow(s,2)*((4. - 2.*delta)*pow(mrho,2) + (2. - 1.*delta)*s) +
                pow(mpion,2)*s*((-4. + 2.*delta)*pow(mrho,2) + (-6. + 3.*delta)*s)) +
             eta1*((2. - 1.*delta)*pow(mpion,6) + (2. - 1.*delta)*pow(mrho,4)*s + (-2. + 1.*delta)*pow(s,3) +
                pow(mpion,4)*((4. - 2.*delta)*pow(mrho,2) + (-6. + 3.*delta)*s) +
                pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) + (-4. + 2.*delta)*pow(mrho,2)*s + (6. - 3.*delta)*pow(s,2))))*
           log(abs(-2.*pow(mpion,2) + s + tmin)))/(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) - 2.*pow(ma1,2)*s + pow(s,2)) +
        (0.25*(0. + 8.*pow(mpion,2)*pow(mrho,6) + 4.*pow(mrho,8) - 8.*pow(mrho,6)*s +
             delta*pow(mrho,4)*(8.*pow(mpion,4) - 8.*pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 16.*s) + 8.*pow(mrho,2)*s +
                8.*pow(s,2)) + pow(delta,2)*pow(mrho,4)*
              (-4.*pow(mpion,4) + 3.*pow(mrho,4) - 2.*pow(mrho,2)*s - 4.*pow(s,2) + pow(mpion,2)*(-6.*pow(mrho,2) + 8.*s)))*
           log(abs(-2.*pow(mpion,2) + 1.*s + 1.*tmin)))/pow(mrho,6) +
        (0.5*(0. + pow(mpion,2)*(4.*pow(mrho,6) - 8.*C4*pow(mrho,8)) - 4.*pow(mrho,6)*s + 8.*C4*pow(mrho,8)*s +
             pow(delta,2)*pow(mrho,4)*(2.*pow(mpion,4) - 1.*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) + 2.*pow(s,2)) +
             delta*pow(mrho,4)*(-4.*pow(mpion,4) + 2.*pow(mrho,2)*s - 4.*pow(s,2) +
                pow(mpion,2)*(-10.*pow(mrho,2) + 4.*C4*pow(mrho,4) + 8.*s) + pow(mrho,4)*(2. - 4.*C4*s)))*
           log(abs(-2.*pow(mpion,2) + 1.*s + 1.*tmin)))/pow(mrho,6)))/
    (16.*Pi*(pow(mpion,4) + pow(pow(mrho,2) - s,2) - 2*pow(mpion,2)*(pow(mrho,2) + s)));


  return sigma / spin_deg_factor;
}

//C13
double total_xsection_C13(double s) {
  double sigma;
  double spin_deg_factor = 3.0;
  double tmin = min_mandelstam_t(s, mpion, mrho, mpion);
  double tmax = max_mandelstam_t(s, mpion, mrho, mpion);

  sigma = (pow(Const,2)*pow(ghat,4)*(0. + (0.03125*pow(eta1 - 1.*eta2,2)*
           (pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) +
                1.*pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
                pow(ma1,2)*pow(mpion,2)*(-2.*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) +
                   2.*pow(mrho,2)*s) + pow(ma1,4)*
                 (4.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-4.*pow(mrho,2) - 4.*s) -
                   2.*pow(mrho,2)*s + 2.*pow(s,2))) +
             eta1*eta2*(-2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) +
                pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) +
                pow(ma1,2)*pow(mpion,2)*(-4.*pow(mrho,4) - 4.*pow(mrho,2)*s +
                   pow(mpion,2)*(4.*pow(mrho,2) + 4.*s)) +
                pow(ma1,4)*(-8.*pow(mpion,4) + 2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 4.*pow(s,2) +
                   pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s))) +
             pow(eta1,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) -
                2.*pow(mpion,2)*pow(mrho,4)*s + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                pow(mpion,4)*(3.*pow(mrho,4) + 2.*pow(mrho,2)*s) +
                pow(ma1,4)*(4.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) -
                   4.*pow(mrho,2)*s + 2.*pow(s,2)) +
                pow(ma1,2)*(pow(mpion,4)*(-6.*pow(mrho,2) - 2.*s) + pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s +
                   pow(mpion,2)*(-4.*pow(mrho,4) + 6.*pow(mrho,2)*s)))))/(1.*pow(ma1,2) - 1.*tmax) +
        (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/
         (1.*pow(mpion,2) - 1.*tmax) - (0.25*pow(-2. + delta,2)*pow(mpion,2)*tmax)/pow(mrho,2) +
        0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-1.*pow(ma1,2) + pow(mrho,2) - 2.*s) +
           eta1*(pow(ma1,2) - 1.*pow(mpion,2) - 2.*pow(mrho,2) + s))*tmax +
        0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta1,2)*
            (3.*pow(ma1,4) + 4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) -
              4.*pow(mrho,2)*s + 2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
           pow(eta2,2)*(3.*pow(ma1,4) + 4.*pow(mpion,4) + pow(mrho,4) +
              pow(mpion,2)*(-4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2) +
              pow(ma1,2)*(-8.*pow(mpion,2) + 4.*pow(mrho,2) + 4.*s)) +
           eta1*eta2*(-6.*pow(ma1,4) - 8.*pow(mpion,4) + 2.*pow(mrho,4) +
              pow(ma1,2)*(16.*pow(mpion,2) - 8.*s) + 4.*pow(mrho,2)*s - 4.*pow(s,2) +
              pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s)))*tmax +
        (2.*(0. - 0.25*pow(mrho,4) - 1.*C4*pow(mrho,6) +
             pow(mpion,2)*(0.75*pow(mrho,2) - 1.*C4*pow(mrho,4)) + 2.*C4*pow(mrho,4)*s +
             pow(delta,2)*(-0.125*pow(mpion,4) - 0.1875*pow(mrho,4) +
                pow(mpion,2)*(0.0625*pow(mrho,2) + 0.0625*s) + 0.1875*pow(mrho,2)*s) +
             delta*(0.25*pow(mpion,4) + 0.5*C4*pow(mrho,6) +
                pow(mpion,2)*(-0.5*pow(mrho,2) + 0.5*C4*pow(mrho,4) - 0.125*s) - 0.375*pow(mrho,2)*s +
                pow(mrho,4)*(0.5 - 1.*C4*s)))*tmax)/pow(mrho,4) -
        (0.0625*(0. + 8.*pow(mpion,2)*pow(mrho,4) - 12.*pow(mrho,6) + 4.*pow(mrho,4)*s +
             delta*pow(mrho,2)*(-16.*pow(mpion,4) - 16.*pow(mpion,2)*pow(mrho,2) - 4.*pow(mrho,4) +
                16.*pow(mrho,2)*s + 4.*pow(s,2)) +
             pow(delta,2)*(8.*pow(mpion,6) + 9.*pow(mrho,6) + pow(mpion,4)*(4.*pow(mrho,2) - 4.*s) -
                13.*pow(mrho,4)*s - 5.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) +
                pow(mpion,2)*(-2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 2.*pow(s,2))))*tmax)/pow(mrho,6) -
        (0.0625*(1.*eta1 - 1.*eta2)*(pow(mrho,2)*
              (eta1*(2.*pow(ma1,2) + 2.*pow(mrho,2)) +
                eta2*(-2.*pow(ma1,2) + 2.*pow(mpion,2) - 8.*pow(mrho,2) + 6.*s)) +
             delta*(eta1*(1.*pow(ma1,4) - 2.*pow(mpion,4) - 3.*pow(mrho,4) +
                   pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 5.000000000000001*pow(s,2) +
                   pow(ma1,2)*(-2.*pow(mpion,2) + 1.*s)) +
                eta2*(-1.*pow(ma1,4) - 4.*pow(mpion,4) + 4.*pow(mrho,4) +
                   pow(mpion,2)*(-1.*pow(mrho,2) - 2.*s) + 1.*pow(mrho,2)*s - 1.*pow(s,2) +
                   pow(ma1,2)*(3.*pow(mpion,2) - 3.*pow(mrho,2) + 2.*s))))*tmax)/pow(mrho,2) -
        (0.5*(pow(mrho,6)*(-1.5 + C4*(-12.*pow(mpion,2) + 6.*s) +
                pow(C4,2)*(-16.*pow(mpion,4) + 16.*pow(mpion,2)*s - 4.*pow(s,2))) +
             pow(delta,2)*(1.*pow(mpion,6) - 2.*pow(mpion,4)*pow(mrho,2) + 0.125*pow(mrho,6) +
                0.25*pow(mrho,4)*s - 0.875*pow(mrho,2)*pow(s,2) + 0.25*pow(s,3) +
                pow(mpion,2)*(-2.5*pow(mrho,4) + 2.25*pow(mrho,2)*s - 0.75*pow(s,2))) +
             delta*pow(mrho,2)*(pow(mpion,4)*(1. + 8.*C4*pow(mrho,2)) + 0.5*pow(s,2) +
                pow(mrho,4)*(1.5 - 5.*C4*s) + pow(mrho,2)*s*(-0.5 + 1.*C4*s) +
                pow(mpion,2)*(6.*C4*pow(mrho,4) - 1.5*s + pow(mrho,2)*(3. - 6.*C4*s))))*tmax)/pow(mrho,6) -
        (0.5*(0. - 4.*C4*pow(mrho,8) - 0.5*pow(mrho,4)*s + pow(mrho,6)*(2. + 2.*C4*s) +
             pow(delta,2)*(-2.*pow(mpion,6) - 2.*pow(mrho,6) + 0.5*pow(mpion,4)*s + 2.125*pow(mrho,4)*s +
                1.25*pow(mrho,2)*pow(s,2) - 0.375*pow(s,3) +
                pow(mpion,2)*(1.5*pow(mrho,4) - 1.5*pow(mrho,2)*s + 1.*pow(s,2))) +
             delta*pow(mrho,2)*(2.*pow(mpion,4) + 2.*C4*pow(mrho,6) - 1.*pow(s,2) +
                pow(mrho,2)*s*(-3. + 1.*C4*s) + pow(mrho,4)*(1. + 1.*C4*s) +
                pow(mpion,2)*(1.*s + pow(mrho,2)*(1. - 2.*C4*s))))*tmax)/pow(mrho,6) +
        (0.0625*(1.*eta1 - 1.*eta2)*(delta*(eta1*
                 (3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(18.*pow(mrho,2) - 12.*s) -
                   8.*pow(mrho,2)*s + 7.*pow(s,2) + pow(ma1,2)*(-10.*pow(mpion,2) - 4.*pow(mrho,2) + 5.*s))\
                 + eta2*(-3.*pow(ma1,4) - 12.*pow(mpion,4) + 2.*pow(mrho,4) +
                   pow(ma1,2)*(11.*pow(mpion,2) - 3.*pow(mrho,2) - 2.*s) + 5.*pow(mrho,2)*s - 3.*pow(s,2) +
                   pow(mpion,2)*(-1.*pow(mrho,2) + 6.*s))) +
             pow(mrho,2)*(eta1*(-8.*C4*pow(ma1,4) - 32.*C4*pow(mpion,4) - 6.*pow(mrho,2) +
                   pow(ma1,2)*(6. + C4*(32.*pow(mpion,2) + 8.*pow(mrho,2) - 16.*s)) + 4.*s +
                   16.*C4*pow(mrho,2)*s - 8.*C4*pow(s,2) + pow(mpion,2)*(-12. - 32.*C4*pow(mrho,2) + 32.*C4*s)
                   ) + eta2*(8.*C4*pow(ma1,4) + 32.*C4*pow(mpion,4) - 4.*pow(mrho,2) - 2.*s +
                   8.*C4*pow(s,2) + pow(mpion,2)*(10. - 16.*C4*pow(mrho,2) - 32.*C4*s) +
                   pow(ma1,2)*(-6. + C4*(-32.*pow(mpion,2) + 8.*pow(mrho,2) + 16.*s)))))*tmax)/pow(mrho,2) +
        0.0625*(-2. + delta)*pow(eta1 - 1.*eta2,2)*pow(tmax,2) +
        (0.1875*(1.3333333333333333*pow(mrho,2) + 5.333333333333333*C4*pow(mrho,4) +
             pow(delta,2)*(1.*pow(mpion,2) + 1.3333333333333333*pow(mrho,2) - 0.3333333333333333*s) +
             delta*(-2.*pow(mpion,2) - 3.3333333333333335*pow(mrho,2) - 2.6666666666666665*C4*pow(mrho,4) +
                0.6666666666666666*s))*pow(tmax,2))/pow(mrho,4) +
        0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) + 2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
           eta1*(pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*pow(tmax,2) -
        (0.375*(0.3333333333333333*pow(mrho,4) - 1.3333333333333333*C4*pow(mrho,6) +
             delta*pow(mrho,2)*(1.3333333333333333*pow(mrho,2) - 0.6666666666666666*C4*pow(mrho,4) +
                pow(mpion,2)*(-0.6666666666666666 + 1.3333333333333333*C4*pow(mrho,2)) - 0.6666666666666666*s) +
             pow(delta,2)*(1.*pow(mpion,4) + 0.25*pow(mrho,4) +
                pow(mpion,2)*(-0.3333333333333333*pow(mrho,2) + 0.6666666666666666*s) -
                0.5833333333333334*pow(s,2)))*pow(tmax,2))/pow(mrho,6) -
        (0.03125*(1.*eta1 - 1.*eta2)*((2.*eta1 - 2.*eta2)*pow(mrho,2) +
             delta*(eta1*(1.*pow(ma1,2) - 2.*pow(mpion,2) + 1.*s) +
                eta2*(-1.*pow(ma1,2) + 3.*pow(mpion,2) - 3.*pow(mrho,2) + 2.*s)))*pow(tmax,2))/pow(mrho,2)\
         + (0.03125*(0. - 4.*pow(mrho,4) + delta*(16.*pow(mrho,4) - 8.*pow(mrho,2)*s) +
             pow(delta,2)*(4.*pow(mpion,4) - 3.*pow(mrho,4) + 2.*pow(mrho,2)*s - 3.*pow(s,2) +
                pow(mpion,2)*(-4.*pow(mrho,2) + 4.*s)))*pow(tmax,2))/pow(mrho,6) +
        (0.25*(C4*pow(mrho,6)*(-6. - 16.*C4*pow(mpion,2) + 8.*C4*s) +
             pow(delta,2)*(1.*pow(mpion,4) - 0.25*pow(mrho,4) + pow(mpion,2)*(-1.75*pow(mrho,2) + 0.5*s) +
                0.5*pow(mrho,2)*s - 0.5*pow(s,2)) +
             delta*pow(mrho,2)*(1.*C4*pow(mrho,4) + pow(mpion,2)*(0.5 + 10.*C4*pow(mrho,2)) - 0.5*s +
                pow(mrho,2)*(2.5 - 4.*C4*s)))*pow(tmax,2))/pow(mrho,6) +
        (0.09375*(1.*eta1 - 1.*eta2)*(delta*(eta2*
                 (-1.*pow(ma1,2) + 3.6666666666666665*pow(mpion,2) - 1.*pow(mrho,2) - 0.6666666666666666*s) +
                eta1*(1.*pow(ma1,2) - 3.3333333333333335*pow(mpion,2) - 1.3333333333333333*pow(mrho,2) +
                   1.6666666666666667*s)) + pow(mrho,2)*
              (eta1*(2. + C4*(-2.6666666666666665*pow(ma1,2) + 10.666666666666666*pow(mpion,2) +
                      2.6666666666666665*pow(mrho,2) - 5.333333333333333*s)) +
                eta2*(-2. + C4*(2.6666666666666665*pow(ma1,2) - 10.666666666666666*pow(mpion,2) +
                      2.6666666666666665*pow(mrho,2) + 5.333333333333333*s))))*pow(tmax,2))/pow(mrho,2) +
        0.010416666666666666*pow(eta1 - 1.*eta2,4)*pow(tmax,3) -
        (0.041666666666666664*delta*(-2. + 1.*delta)*pow(tmax,3))/pow(mrho,4) -
        (0.020833333333333332*delta*pow(1.*eta1 - 1.*eta2,2)*pow(tmax,3))/pow(mrho,2) -
        (0.16666666666666666*pow(1.*eta1 - 1.*eta2,2)*(-0.375*delta + 1.*C4*pow(mrho,2))*pow(tmax,3))/
         pow(mrho,2) + (0.10416666666666666*delta*
           (-0.8*pow(mrho,2) + delta*(0.4*pow(mpion,2) + 1.*pow(mrho,2) - 0.6*s))*pow(tmax,3))/pow(mrho,6)\
         + (0.16666666666666666*delta*(0. - 0.75*delta*pow(mrho,2) + 1.*C4*pow(mrho,4) + 0.625*delta*s)*
           pow(tmax,3))/pow(mrho,6) - (0.041666666666666664*
           (12.*C4*delta*pow(mrho,4) - 16.*pow(C4,2)*pow(mrho,6) +
             pow(delta,2)*(1.*pow(mpion,2) - 2.5*pow(mrho,2) + 1.*s))*pow(tmax,3))/pow(mrho,6) +
        (0. - 0.5000000000000001*pow(2. - 1.*delta,2)*pow(mpion,4)*pow(mrho,6) +
           0.25*pow(2. - 1.*delta,2)*pow(mrho,10) -
           0.5000000000000001*pow(2. - 1.*delta,2)*pow(mrho,6)*pow(s,2) +
           pow(2. - 1.*delta,2)*pow(mpion,2)*pow(mrho,6)*(-1.*pow(mrho,2) + 1.*s))/
         (pow(mrho,6)*(-2.*pow(mpion,2) + 1.*s + 1.*tmax)) -
        (0.0625*(1.*eta1 - 1.*eta2)*(2.*pow(mrho,2) +
             delta*(1.*pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s))*
           (eta2*(-1.*pow(ma1,6) + pow(mpion,2)*(4.*pow(mpion,2) - 1.*s)*s +
                pow(ma1,4)*(3.*pow(mpion,2) - 4.*pow(mrho,2) + 2.*s) +
                pow(ma1,2)*(-4.*pow(mpion,4) - 2.*pow(mpion,2)*s + (4.*pow(mrho,2) - 1.*s)*s)) +
             eta1*(1.*pow(ma1,6) + 8.*pow(mpion,6) + pow(mpion,4)*(-4.*pow(mrho,2) - 6.*s) +
                pow(mpion,2)*(4.*pow(mrho,2) - 2.*s)*s +
                pow(ma1,4)*(-2.*pow(mpion,2) + 1.*pow(mrho,2) + 1.*s) +
                s*(2.*pow(mrho,4) - 3.*pow(mrho,2)*s + 1.*pow(s,2)) +
                pow(ma1,2)*(-2.*pow(mpion,4) - 2.*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) -
                   2.*pow(mrho,2)*s + 5.*pow(s,2))))*log(abs(-pow(ma1,2) + tmax)))/
         (pow(mrho,2)*(pow(ma1,2) - 2.*pow(mpion,2) + s)) +
        (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*
              (-1.*pow(ma1,6) + pow(mpion,6) - 1.*pow(mpion,4)*pow(mrho,2) +
                pow(ma1,4)*(pow(mpion,2) + pow(mrho,2) - 2.*s) +
                pow(ma1,2)*(3.*pow(mpion,4) - 2.*pow(mpion,2)*s)) +
             eta1*(pow(ma1,6) + pow(ma1,4)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + s) +
                pow(mpion,2)*(-4.*pow(mpion,4) - 1.*pow(mrho,4) - 1.*pow(mrho,2)*s +
                   pow(mpion,2)*(3.*pow(mrho,2) + s)) +
                pow(ma1,2)*(pow(mpion,4) + pow(mrho,4) - 1.*pow(mrho,2)*s +
                   pow(mpion,2)*(pow(mrho,2) + 2.*s))))*log(abs(-pow(ma1,2) + tmax)))/
         (pow(ma1,2) - 1.*pow(mpion,2)) + 0.0625*pow(eta1 - 1.*eta2,2)*
         (pow(eta1,2)*(2.*pow(ma1,6) + pow(mpion,4)*(-3.*pow(mrho,2) - 1.*s) +
              pow(mrho,2)*(pow(mrho,2) - 1.*s)*s + pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
              pow(mpion,2)*(-2.*pow(mrho,4) + 3.*pow(mrho,2)*s) +
              pow(ma1,2)*(4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) -
                 4.*pow(mrho,2)*s + 2.*pow(s,2))) +
           pow(eta2,2)*(2.*pow(ma1,6) + pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
              pow(mpion,2)*(-1.*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 1.*s) + pow(mrho,2)*s) +
              pow(ma1,2)*(4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-4.*pow(mrho,2) - 4.*s) -
                 2.*pow(mrho,2)*s + 2.*pow(s,2))) +
           eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) +
              pow(mpion,2)*(-2.*pow(mrho,4) - 2.*pow(mrho,2)*s + pow(mpion,2)*(2.*pow(mrho,2) + 2.*s)) +
              pow(ma1,2)*(-8.*pow(mpion,4) + 2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 4.*pow(s,2) +
                 pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s))))*log(abs(-pow(ma1,2) + tmax)) +
        (0.0625*(1.*eta1 - 1.*eta2)*(delta*(eta1*
                 (3.*pow(ma1,6) + 8.*pow(mpion,6) + pow(mpion,4)*(-12.*pow(mrho,2) - 6.*s) +
                   pow(ma1,4)*(-10.*pow(mpion,2) - 4.*pow(mrho,2) + 5.*s) +
                   pow(mpion,2)*(-4.*pow(mrho,4) + 10.*pow(mrho,2)*s - 2.*pow(s,2)) +
                   s*(3.*pow(mrho,4) - 4.*pow(mrho,2)*s + pow(s,2)) +
                   pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(18.*pow(mrho,2) - 12.*s) -
                      8.*pow(mrho,2)*s + 7.*pow(s,2))) +
                eta2*(-3.*pow(ma1,6) + pow(ma1,4)*(11.*pow(mpion,2) - 3.*pow(mrho,2) - 2.*s) +
                   pow(mpion,2)*(-2.*pow(mrho,4) + pow(mrho,2)*s - 1.*pow(s,2) +
                      pow(mpion,2)*(-2.*pow(mrho,2) + 4.*s)) +
                   pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 5.*pow(mrho,2)*s - 3.*pow(s,2) +
                      pow(mpion,2)*(-1.*pow(mrho,2) + 6.*s)))) +
             pow(mrho,2)*(eta2*(8.*C4*pow(ma1,6) +
                   pow(mpion,2)*((4. + 8.*C4*pow(mpion,2))*pow(mrho,2) - 2.*s) +
                   pow(ma1,4)*(-6. + C4*(-32.*pow(mpion,2) + 8.*pow(mrho,2) + 16.*s)) +
                   pow(ma1,2)*(32.*C4*pow(mpion,4) - 4.*pow(mrho,2) +
                      pow(mpion,2)*(10. - 16.*C4*pow(mrho,2) - 32.*C4*s) + s*(-2. + 8.*C4*s))) +
                eta1*(-8.*C4*pow(ma1,6) + pow(mpion,4)*(4. + 24.*C4*pow(mrho,2)) +
                   pow(ma1,4)*(6. + C4*(32.*pow(mpion,2) + 8.*pow(mrho,2) - 16.*s)) +
                   s*(-2.*pow(mrho,2) + 2.*s) + pow(mpion,2)*(-4.*s + pow(mrho,2)*(8. - 16.*C4*s)) +
                   pow(ma1,2)*(-32.*C4*pow(mpion,4) + s*(4. - 8.*C4*s) + pow(mrho,2)*(-6. + 16.*C4*s) +
                      pow(mpion,2)*(-12. - 32.*C4*pow(mrho,2) + 32.*C4*s)))))*log(abs(-pow(ma1,2) + tmax)))/
         pow(mrho,2) + 0.5*pow(-2. + delta,2)*pow(mpion,2)*log(abs(-pow(mpion,2) + tmax)) -
        (0.25*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(2.*pow(mpion,6) - 2.*pow(mpion,4)*s) +
             eta1*(-2.*pow(mpion,6) - 1.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(pow(mrho,2) + 2.*s)))*
           log(abs(-pow(mpion,2) + tmax)))/(pow(ma1,2) - 1.*pow(mpion,2)) -
        (0.125*(0. - 32.*C4*pow(mpion,6)*pow(mrho,4) - 8.*pow(mrho,8) + 8.*pow(mrho,6)*s +
             pow(mpion,4)*pow(mrho,4)*(16. + 64.*C4*s) +
             pow(mpion,2)*pow(mrho,4)*(24.*pow(mrho,2) + s*(-16. - 32.*C4*s)) +
             pow(delta,2)*pow(mrho,2)*(-4.*pow(mpion,6) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s +
                pow(mpion,4)*(4.*pow(mrho,2) + 8.*s) +
                pow(mpion,2)*(6.*pow(mrho,4) - 4.*pow(mrho,2)*s - 4.000000000000001*pow(s,2))) +
             delta*pow(mrho,2)*(8.*pow(mrho,6) + pow(mpion,6)*(8. + 16.*C4*pow(mrho,2)) -
                8.*pow(mrho,4)*s + pow(mpion,4)*(-16.*s + pow(mrho,2)*(-15.999999999999996 - 32.*C4*s)) +
                pow(mpion,2)*(-24.*pow(mrho,4) + 8.*pow(s,2) + pow(mrho,2)*s*(16. + 16.*C4*s))))*
           log(abs(-pow(mpion,2) + tmax)))/(pow(mrho,4)*(-1.*pow(mpion,2) + 1.*s)) -
        (0.25*(1.*eta1 - 1.*eta2)*(eta2*((2. - 1.*delta)*pow(mpion,6) +
                pow(mpion,2)*s*((-12. + 6.*delta)*pow(mrho,2) + (6. - 3.*delta)*s) +
                pow(s,2)*((4. - 2.*delta)*pow(mrho,2) + (-2. + 1.*delta)*s) +
                pow(mpion,4)*((8. - 4.*delta)*pow(mrho,2) + (-6. + 3.*delta)*s)) +
             eta1*((-2. + 1.*delta)*pow(mpion,6) + (-2. + 1.*delta)*pow(mrho,4)*s + (2. - 1.*delta)*pow(s,3) +
                pow(mpion,4)*((-4. + 2.*delta)*pow(mrho,2) + (6. - 3.*delta)*s) +
                pow(mpion,2)*((2. - 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s +
                   (-6. + 3.*delta)*pow(s,2))))*log(abs(-2.*pow(mpion,2) + s + tmax)))/
         (pow(ma1,2) - 2.*pow(mpion,2) + s) +
        (0.125*(0. + (32. - 31.999999999999993*delta + 8.*pow(delta,2))*pow(mpion,4)*pow(mrho,4) -
             2.0000000000000004*pow(2. - 1.*delta,2)*pow(mrho,8) +
             pow(mpion,2)*pow(mrho,4)*(8.000000000000002*pow(2. - 1.*delta,2)*pow(mrho,2) +
                (-32. + 31.999999999999996*delta - 8.*pow(delta,2))*s))*log(abs(-2.*pow(mpion,2) + s + tmax)))/
         (pow(mrho,4)*(-1.*pow(mpion,2) + 1.*s)) +
        (0.25*(0. + 8.*pow(mpion,2)*pow(mrho,6) + 4.*pow(mrho,8) - 8.*pow(mrho,6)*s +
             delta*pow(mrho,4)*(8.*pow(mpion,4) - 8.*pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 16.*s) +
                8.*pow(mrho,2)*s + 8.*pow(s,2)) +
             pow(delta,2)*pow(mrho,4)*(-4.*pow(mpion,4) + 3.*pow(mrho,4) - 2.*pow(mrho,2)*s -
                4.*pow(s,2) + pow(mpion,2)*(-6.*pow(mrho,2) + 8.*s)))*log(abs(-2.*pow(mpion,2) + s + tmax)))/
         pow(mrho,6) - (0.5*(0. + pow(mpion,2)*(4.*pow(mrho,6) - 8.*C4*pow(mrho,8)) - 4.*pow(mrho,6)*s +
             8.*C4*pow(mrho,8)*s + pow(delta,2)*pow(mrho,4)*
              (-2.*pow(mpion,4) + 1.*pow(mrho,4) - 2.*pow(s,2) + pow(mpion,2)*(-4.*pow(mrho,2) + 4.*s)) +
             delta*pow(mrho,4)*(4.*pow(mpion,4) +
                pow(mpion,2)*(6.*pow(mrho,2) + 4.*C4*pow(mrho,4) - 8.*s) + 2.*pow(mrho,2)*s + 4.*pow(s,2) +
                pow(mrho,4)*(-2. - 4.*C4*s)))*log(abs(-2.*pow(mpion,2) + s + tmax)))/pow(mrho,6)))/
    (16.*Pi*(pow(mpion,4) + pow(pow(mrho,2) - s,2) - 2*pow(mpion,2)*(pow(mrho,2) + s))) -
   (pow(Const,2)*pow(ghat,4)*(0. + (0.03125*pow(eta1 - 1.*eta2,2)*
           (pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) +
                1.*pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
                pow(ma1,2)*pow(mpion,2)*(-2.*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) +
                   2.*pow(mrho,2)*s) + pow(ma1,4)*
                 (4.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-4.*pow(mrho,2) - 4.*s) -
                   2.*pow(mrho,2)*s + 2.*pow(s,2))) +
             eta1*eta2*(-2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) +
                pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) +
                pow(ma1,2)*pow(mpion,2)*(-4.*pow(mrho,4) - 4.*pow(mrho,2)*s +
                   pow(mpion,2)*(4.*pow(mrho,2) + 4.*s)) +
                pow(ma1,4)*(-8.*pow(mpion,4) + 2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 4.*pow(s,2) +
                   pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s))) +
             pow(eta1,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) -
                2.*pow(mpion,2)*pow(mrho,4)*s + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
                pow(mpion,4)*(3.*pow(mrho,4) + 2.*pow(mrho,2)*s) +
                pow(ma1,4)*(4.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) -
                   4.*pow(mrho,2)*s + 2.*pow(s,2)) +
                pow(ma1,2)*(pow(mpion,4)*(-6.*pow(mrho,2) - 2.*s) + pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s +
                   pow(mpion,2)*(-4.*pow(mrho,4) + 6.*pow(mrho,2)*s)))))/(1.*pow(ma1,2) - 1.*tmin) +
        (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/
         (1.*pow(mpion,2) - 1.*tmin) - (0.25*pow(-2. + delta,2)*pow(mpion,2)*tmin)/pow(mrho,2) +
        0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-1.*pow(ma1,2) + pow(mrho,2) - 2.*s) +
           eta1*(pow(ma1,2) - 1.*pow(mpion,2) - 2.*pow(mrho,2) + s))*tmin +
        0.03125*pow(eta1 - 1.*eta2,2)*(pow(eta1,2)*
            (3.*pow(ma1,4) + 4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) -
              4.*pow(mrho,2)*s + 2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
           pow(eta2,2)*(3.*pow(ma1,4) + 4.*pow(mpion,4) + pow(mrho,4) +
              pow(mpion,2)*(-4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2) +
              pow(ma1,2)*(-8.*pow(mpion,2) + 4.*pow(mrho,2) + 4.*s)) +
           eta1*eta2*(-6.*pow(ma1,4) - 8.*pow(mpion,4) + 2.*pow(mrho,4) +
              pow(ma1,2)*(16.*pow(mpion,2) - 8.*s) + 4.*pow(mrho,2)*s - 4.*pow(s,2) +
              pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s)))*tmin +
        (2.*(0. - 0.25*pow(mrho,4) - 1.*C4*pow(mrho,6) +
             pow(mpion,2)*(0.75*pow(mrho,2) - 1.*C4*pow(mrho,4)) + 2.*C4*pow(mrho,4)*s +
             pow(delta,2)*(-0.125*pow(mpion,4) - 0.1875*pow(mrho,4) +
                pow(mpion,2)*(0.0625*pow(mrho,2) + 0.0625*s) + 0.1875*pow(mrho,2)*s) +
             delta*(0.25*pow(mpion,4) + 0.5*C4*pow(mrho,6) +
                pow(mpion,2)*(-0.5*pow(mrho,2) + 0.5*C4*pow(mrho,4) - 0.125*s) - 0.375*pow(mrho,2)*s +
                pow(mrho,4)*(0.5 - 1.*C4*s)))*tmin)/pow(mrho,4) -
        (0.0625*(0. + 8.*pow(mpion,2)*pow(mrho,4) - 12.*pow(mrho,6) + 4.*pow(mrho,4)*s +
             delta*pow(mrho,2)*(-16.*pow(mpion,4) - 16.*pow(mpion,2)*pow(mrho,2) - 4.*pow(mrho,4) +
                16.*pow(mrho,2)*s + 4.*pow(s,2)) +
             pow(delta,2)*(8.*pow(mpion,6) + 9.*pow(mrho,6) + pow(mpion,4)*(4.*pow(mrho,2) - 4.*s) -
                13.*pow(mrho,4)*s - 5.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) +
                pow(mpion,2)*(-2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 2.*pow(s,2))))*tmin)/pow(mrho,6) -
        (0.0625*(1.*eta1 - 1.*eta2)*(pow(mrho,2)*
              (eta1*(2.*pow(ma1,2) + 2.*pow(mrho,2)) +
                eta2*(-2.*pow(ma1,2) + 2.*pow(mpion,2) - 8.*pow(mrho,2) + 6.*s)) +
             delta*(eta1*(1.*pow(ma1,4) - 2.*pow(mpion,4) - 3.*pow(mrho,4) +
                   pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 5.000000000000001*pow(s,2) +
                   pow(ma1,2)*(-2.*pow(mpion,2) + 1.*s)) +
                eta2*(-1.*pow(ma1,4) - 4.*pow(mpion,4) + 4.*pow(mrho,4) +
                   pow(mpion,2)*(-1.*pow(mrho,2) - 2.*s) + 1.*pow(mrho,2)*s - 1.*pow(s,2) +
                   pow(ma1,2)*(3.*pow(mpion,2) - 3.*pow(mrho,2) + 2.*s))))*tmin)/pow(mrho,2) -
        (0.5*(pow(mrho,6)*(-1.5 + C4*(-12.*pow(mpion,2) + 6.*s) +
                pow(C4,2)*(-16.*pow(mpion,4) + 16.*pow(mpion,2)*s - 4.*pow(s,2))) +
             pow(delta,2)*(1.*pow(mpion,6) - 2.*pow(mpion,4)*pow(mrho,2) + 0.125*pow(mrho,6) +
                0.25*pow(mrho,4)*s - 0.875*pow(mrho,2)*pow(s,2) + 0.25*pow(s,3) +
                pow(mpion,2)*(-2.5*pow(mrho,4) + 2.25*pow(mrho,2)*s - 0.75*pow(s,2))) +
             delta*pow(mrho,2)*(pow(mpion,4)*(1. + 8.*C4*pow(mrho,2)) + 0.5*pow(s,2) +
                pow(mrho,4)*(1.5 - 5.*C4*s) + pow(mrho,2)*s*(-0.5 + 1.*C4*s) +
                pow(mpion,2)*(6.*C4*pow(mrho,4) - 1.5*s + pow(mrho,2)*(3. - 6.*C4*s))))*tmin)/pow(mrho,6) -
        (0.5*(0. - 4.*C4*pow(mrho,8) - 0.5*pow(mrho,4)*s + pow(mrho,6)*(2. + 2.*C4*s) +
             pow(delta,2)*(-2.*pow(mpion,6) - 2.*pow(mrho,6) + 0.5*pow(mpion,4)*s + 2.125*pow(mrho,4)*s +
                1.25*pow(mrho,2)*pow(s,2) - 0.375*pow(s,3) +
                pow(mpion,2)*(1.5*pow(mrho,4) - 1.5*pow(mrho,2)*s + 1.*pow(s,2))) +
             delta*pow(mrho,2)*(2.*pow(mpion,4) + 2.*C4*pow(mrho,6) - 1.*pow(s,2) +
                pow(mrho,2)*s*(-3. + 1.*C4*s) + pow(mrho,4)*(1. + 1.*C4*s) +
                pow(mpion,2)*(1.*s + pow(mrho,2)*(1. - 2.*C4*s))))*tmin)/pow(mrho,6) +
        (0.0625*(1.*eta1 - 1.*eta2)*(delta*(eta1*
                 (3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(18.*pow(mrho,2) - 12.*s) -
                   8.*pow(mrho,2)*s + 7.*pow(s,2) + pow(ma1,2)*(-10.*pow(mpion,2) - 4.*pow(mrho,2) + 5.*s))\
                 + eta2*(-3.*pow(ma1,4) - 12.*pow(mpion,4) + 2.*pow(mrho,4) +
                   pow(ma1,2)*(11.*pow(mpion,2) - 3.*pow(mrho,2) - 2.*s) + 5.*pow(mrho,2)*s - 3.*pow(s,2) +
                   pow(mpion,2)*(-1.*pow(mrho,2) + 6.*s))) +
             pow(mrho,2)*(eta1*(-8.*C4*pow(ma1,4) - 32.*C4*pow(mpion,4) - 6.*pow(mrho,2) +
                   pow(ma1,2)*(6. + C4*(32.*pow(mpion,2) + 8.*pow(mrho,2) - 16.*s)) + 4.*s +
                   16.*C4*pow(mrho,2)*s - 8.*C4*pow(s,2) + pow(mpion,2)*(-12. - 32.*C4*pow(mrho,2) + 32.*C4*s)
                   ) + eta2*(8.*C4*pow(ma1,4) + 32.*C4*pow(mpion,4) - 4.*pow(mrho,2) - 2.*s +
                   8.*C4*pow(s,2) + pow(mpion,2)*(10. - 16.*C4*pow(mrho,2) - 32.*C4*s) +
                   pow(ma1,2)*(-6. + C4*(-32.*pow(mpion,2) + 8.*pow(mrho,2) + 16.*s)))))*tmin)/pow(mrho,2) +
        0.0625*(-2. + delta)*pow(eta1 - 1.*eta2,2)*pow(tmin,2) +
        (0.1875*(1.3333333333333333*pow(mrho,2) + 5.333333333333333*C4*pow(mrho,4) +
             pow(delta,2)*(1.*pow(mpion,2) + 1.3333333333333333*pow(mrho,2) - 0.3333333333333333*s) +
             delta*(-2.*pow(mpion,2) - 3.3333333333333335*pow(mrho,2) - 2.6666666666666665*C4*pow(mrho,4) +
                0.6666666666666666*s))*pow(tmin,2))/pow(mrho,4) +
        0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) + 2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
           eta1*(pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*pow(tmin,2) -
        (0.375*(0.3333333333333333*pow(mrho,4) - 1.3333333333333333*C4*pow(mrho,6) +
             delta*pow(mrho,2)*(1.3333333333333333*pow(mrho,2) - 0.6666666666666666*C4*pow(mrho,4) +
                pow(mpion,2)*(-0.6666666666666666 + 1.3333333333333333*C4*pow(mrho,2)) - 0.6666666666666666*s) +
             pow(delta,2)*(1.*pow(mpion,4) + 0.25*pow(mrho,4) +
                pow(mpion,2)*(-0.3333333333333333*pow(mrho,2) + 0.6666666666666666*s) -
                0.5833333333333334*pow(s,2)))*pow(tmin,2))/pow(mrho,6) -
        (0.03125*(1.*eta1 - 1.*eta2)*((2.*eta1 - 2.*eta2)*pow(mrho,2) +
             delta*(eta1*(1.*pow(ma1,2) - 2.*pow(mpion,2) + 1.*s) +
                eta2*(-1.*pow(ma1,2) + 3.*pow(mpion,2) - 3.*pow(mrho,2) + 2.*s)))*pow(tmin,2))/pow(mrho,2)\
         + (0.03125*(0. - 4.*pow(mrho,4) + delta*(16.*pow(mrho,4) - 8.*pow(mrho,2)*s) +
             pow(delta,2)*(4.*pow(mpion,4) - 3.*pow(mrho,4) + 2.*pow(mrho,2)*s - 3.*pow(s,2) +
                pow(mpion,2)*(-4.*pow(mrho,2) + 4.*s)))*pow(tmin,2))/pow(mrho,6) +
        (0.25*(C4*pow(mrho,6)*(-6. - 16.*C4*pow(mpion,2) + 8.*C4*s) +
             pow(delta,2)*(1.*pow(mpion,4) - 0.25*pow(mrho,4) + pow(mpion,2)*(-1.75*pow(mrho,2) + 0.5*s) +
                0.5*pow(mrho,2)*s - 0.5*pow(s,2)) +
             delta*pow(mrho,2)*(1.*C4*pow(mrho,4) + pow(mpion,2)*(0.5 + 10.*C4*pow(mrho,2)) - 0.5*s +
                pow(mrho,2)*(2.5 - 4.*C4*s)))*pow(tmin,2))/pow(mrho,6) +
        (0.09375*(1.*eta1 - 1.*eta2)*(delta*(eta2*
                 (-1.*pow(ma1,2) + 3.6666666666666665*pow(mpion,2) - 1.*pow(mrho,2) - 0.6666666666666666*s) +
                eta1*(1.*pow(ma1,2) - 3.3333333333333335*pow(mpion,2) - 1.3333333333333333*pow(mrho,2) +
                   1.6666666666666667*s)) + pow(mrho,2)*
              (eta1*(2. + C4*(-2.6666666666666665*pow(ma1,2) + 10.666666666666666*pow(mpion,2) +
                      2.6666666666666665*pow(mrho,2) - 5.333333333333333*s)) +
                eta2*(-2. + C4*(2.6666666666666665*pow(ma1,2) - 10.666666666666666*pow(mpion,2) +
                      2.6666666666666665*pow(mrho,2) + 5.333333333333333*s))))*pow(tmin,2))/pow(mrho,2) +
        0.010416666666666666*pow(eta1 - 1.*eta2,4)*pow(tmin,3) -
        (0.041666666666666664*delta*(-2. + 1.*delta)*pow(tmin,3))/pow(mrho,4) -
        (0.020833333333333332*delta*pow(1.*eta1 - 1.*eta2,2)*pow(tmin,3))/pow(mrho,2) -
        (0.16666666666666666*pow(1.*eta1 - 1.*eta2,2)*(-0.375*delta + 1.*C4*pow(mrho,2))*pow(tmin,3))/
         pow(mrho,2) + (0.10416666666666666*delta*
           (-0.8*pow(mrho,2) + delta*(0.4*pow(mpion,2) + 1.*pow(mrho,2) - 0.6*s))*pow(tmin,3))/pow(mrho,6)\
         + (0.16666666666666666*delta*(0. - 0.75*delta*pow(mrho,2) + 1.*C4*pow(mrho,4) + 0.625*delta*s)*
           pow(tmin,3))/pow(mrho,6) - (0.041666666666666664*
           (12.*C4*delta*pow(mrho,4) - 16.*pow(C4,2)*pow(mrho,6) +
             pow(delta,2)*(1.*pow(mpion,2) - 2.5*pow(mrho,2) + 1.*s))*pow(tmin,3))/pow(mrho,6) +
        (0. - 0.5000000000000001*pow(2. - 1.*delta,2)*pow(mpion,4)*pow(mrho,6) +
           0.25*pow(2. - 1.*delta,2)*pow(mrho,10) -
           0.5000000000000001*pow(2. - 1.*delta,2)*pow(mrho,6)*pow(s,2) +
           pow(2. - 1.*delta,2)*pow(mpion,2)*pow(mrho,6)*(-1.*pow(mrho,2) + 1.*s))/
         (pow(mrho,6)*(-2.*pow(mpion,2) + 1.*s + 1.*tmin)) -
        (0.0625*(1.*eta1 - 1.*eta2)*(2.*pow(mrho,2) +
             delta*(1.*pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s))*
           (eta2*(-1.*pow(ma1,6) + pow(mpion,2)*(4.*pow(mpion,2) - 1.*s)*s +
                pow(ma1,4)*(3.*pow(mpion,2) - 4.*pow(mrho,2) + 2.*s) +
                pow(ma1,2)*(-4.*pow(mpion,4) - 2.*pow(mpion,2)*s + (4.*pow(mrho,2) - 1.*s)*s)) +
             eta1*(1.*pow(ma1,6) + 8.*pow(mpion,6) + pow(mpion,4)*(-4.*pow(mrho,2) - 6.*s) +
                pow(mpion,2)*(4.*pow(mrho,2) - 2.*s)*s +
                pow(ma1,4)*(-2.*pow(mpion,2) + 1.*pow(mrho,2) + 1.*s) +
                s*(2.*pow(mrho,4) - 3.*pow(mrho,2)*s + 1.*pow(s,2)) +
                pow(ma1,2)*(-2.*pow(mpion,4) - 2.*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) -
                   2.*pow(mrho,2)*s + 5.*pow(s,2))))*log(abs(-pow(ma1,2) + tmin)))/
         (pow(mrho,2)*(pow(ma1,2) - 2.*pow(mpion,2) + s)) +
        (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*
              (-1.*pow(ma1,6) + pow(mpion,6) - 1.*pow(mpion,4)*pow(mrho,2) +
                pow(ma1,4)*(pow(mpion,2) + pow(mrho,2) - 2.*s) +
                pow(ma1,2)*(3.*pow(mpion,4) - 2.*pow(mpion,2)*s)) +
             eta1*(pow(ma1,6) + pow(ma1,4)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + s) +
                pow(mpion,2)*(-4.*pow(mpion,4) - 1.*pow(mrho,4) - 1.*pow(mrho,2)*s +
                   pow(mpion,2)*(3.*pow(mrho,2) + s)) +
                pow(ma1,2)*(pow(mpion,4) + pow(mrho,4) - 1.*pow(mrho,2)*s +
                   pow(mpion,2)*(pow(mrho,2) + 2.*s))))*log(abs(-pow(ma1,2) + tmin)))/
         (pow(ma1,2) - 1.*pow(mpion,2)) + 0.0625*pow(eta1 - 1.*eta2,2)*
         (pow(eta1,2)*(2.*pow(ma1,6) + pow(mpion,4)*(-3.*pow(mrho,2) - 1.*s) +
              pow(mrho,2)*(pow(mrho,2) - 1.*s)*s + pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
              pow(mpion,2)*(-2.*pow(mrho,4) + 3.*pow(mrho,2)*s) +
              pow(ma1,2)*(4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 4.*s) -
                 4.*pow(mrho,2)*s + 2.*pow(s,2))) +
           pow(eta2,2)*(2.*pow(ma1,6) + pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
              pow(mpion,2)*(-1.*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 1.*s) + pow(mrho,2)*s) +
              pow(ma1,2)*(4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-4.*pow(mrho,2) - 4.*s) -
                 2.*pow(mrho,2)*s + 2.*pow(s,2))) +
           eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) +
              pow(mpion,2)*(-2.*pow(mrho,4) - 2.*pow(mrho,2)*s + pow(mpion,2)*(2.*pow(mrho,2) + 2.*s)) +
              pow(ma1,2)*(-8.*pow(mpion,4) + 2.*pow(mrho,4) + 4.*pow(mrho,2)*s - 4.*pow(s,2) +
                 pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s))))*log(abs(-pow(ma1,2) + tmin)) +
        (0.0625*(1.*eta1 - 1.*eta2)*(delta*(eta1*
                 (3.*pow(ma1,6) + 8.*pow(mpion,6) + pow(mpion,4)*(-12.*pow(mrho,2) - 6.*s) +
                   pow(ma1,4)*(-10.*pow(mpion,2) - 4.*pow(mrho,2) + 5.*s) +
                   pow(mpion,2)*(-4.*pow(mrho,4) + 10.*pow(mrho,2)*s - 2.*pow(s,2)) +
                   s*(3.*pow(mrho,4) - 4.*pow(mrho,2)*s + pow(s,2)) +
                   pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(18.*pow(mrho,2) - 12.*s) -
                      8.*pow(mrho,2)*s + 7.*pow(s,2))) +
                eta2*(-3.*pow(ma1,6) + pow(ma1,4)*(11.*pow(mpion,2) - 3.*pow(mrho,2) - 2.*s) +
                   pow(mpion,2)*(-2.*pow(mrho,4) + pow(mrho,2)*s - 1.*pow(s,2) +
                      pow(mpion,2)*(-2.*pow(mrho,2) + 4.*s)) +
                   pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 5.*pow(mrho,2)*s - 3.*pow(s,2) +
                      pow(mpion,2)*(-1.*pow(mrho,2) + 6.*s)))) +
             pow(mrho,2)*(eta2*(8.*C4*pow(ma1,6) +
                   pow(mpion,2)*((4. + 8.*C4*pow(mpion,2))*pow(mrho,2) - 2.*s) +
                   pow(ma1,4)*(-6. + C4*(-32.*pow(mpion,2) + 8.*pow(mrho,2) + 16.*s)) +
                   pow(ma1,2)*(32.*C4*pow(mpion,4) - 4.*pow(mrho,2) +
                      pow(mpion,2)*(10. - 16.*C4*pow(mrho,2) - 32.*C4*s) + s*(-2. + 8.*C4*s))) +
                eta1*(-8.*C4*pow(ma1,6) + pow(mpion,4)*(4. + 24.*C4*pow(mrho,2)) +
                   pow(ma1,4)*(6. + C4*(32.*pow(mpion,2) + 8.*pow(mrho,2) - 16.*s)) +
                   s*(-2.*pow(mrho,2) + 2.*s) + pow(mpion,2)*(-4.*s + pow(mrho,2)*(8. - 16.*C4*s)) +
                   pow(ma1,2)*(-32.*C4*pow(mpion,4) + s*(4. - 8.*C4*s) + pow(mrho,2)*(-6. + 16.*C4*s) +
                      pow(mpion,2)*(-12. - 32.*C4*pow(mrho,2) + 32.*C4*s)))))*log(abs(-pow(ma1,2) + tmin)))/
         pow(mrho,2) + 0.5*pow(-2. + delta,2)*pow(mpion,2)*log(abs(-pow(mpion,2) + tmin)) -
        (0.25*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(2.*pow(mpion,6) - 2.*pow(mpion,4)*s) +
             eta1*(-2.*pow(mpion,6) - 1.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(pow(mrho,2) + 2.*s)))*
           log(abs(-pow(mpion,2) + tmin)))/(pow(ma1,2) - 1.*pow(mpion,2)) -
        (0.125*(0. - 32.*C4*pow(mpion,6)*pow(mrho,4) - 8.*pow(mrho,8) + 8.*pow(mrho,6)*s +
             pow(mpion,4)*pow(mrho,4)*(16. + 64.*C4*s) +
             pow(mpion,2)*pow(mrho,4)*(24.*pow(mrho,2) + s*(-16. - 32.*C4*s)) +
             pow(delta,2)*pow(mrho,2)*(-4.*pow(mpion,6) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s +
                pow(mpion,4)*(4.*pow(mrho,2) + 8.*s) +
                pow(mpion,2)*(6.*pow(mrho,4) - 4.*pow(mrho,2)*s - 4.000000000000001*pow(s,2))) +
             delta*pow(mrho,2)*(8.*pow(mrho,6) + pow(mpion,6)*(8. + 16.*C4*pow(mrho,2)) -
                8.*pow(mrho,4)*s + pow(mpion,4)*(-16.*s + pow(mrho,2)*(-15.999999999999996 - 32.*C4*s)) +
                pow(mpion,2)*(-24.*pow(mrho,4) + 8.*pow(s,2) + pow(mrho,2)*s*(16. + 16.*C4*s))))*
           log(abs(-pow(mpion,2) + tmin)))/(pow(mrho,4)*(-1.*pow(mpion,2) + 1.*s)) -
        (0.25*(1.*eta1 - 1.*eta2)*(eta2*((2. - 1.*delta)*pow(mpion,6) +
                pow(mpion,2)*s*((-12. + 6.*delta)*pow(mrho,2) + (6. - 3.*delta)*s) +
                pow(s,2)*((4. - 2.*delta)*pow(mrho,2) + (-2. + 1.*delta)*s) +
                pow(mpion,4)*((8. - 4.*delta)*pow(mrho,2) + (-6. + 3.*delta)*s)) +
             eta1*((-2. + 1.*delta)*pow(mpion,6) + (-2. + 1.*delta)*pow(mrho,4)*s + (2. - 1.*delta)*pow(s,3) +
                pow(mpion,4)*((-4. + 2.*delta)*pow(mrho,2) + (6. - 3.*delta)*s) +
                pow(mpion,2)*((2. - 1.*delta)*pow(mrho,4) + (4. - 2.*delta)*pow(mrho,2)*s +
                   (-6. + 3.*delta)*pow(s,2))))*log(abs(-2.*pow(mpion,2) + s + tmin)))/
         (pow(ma1,2) - 2.*pow(mpion,2) + s) +
        (0.125*(0. + (32. - 31.999999999999993*delta + 8.*pow(delta,2))*pow(mpion,4)*pow(mrho,4) -
             2.0000000000000004*pow(2. - 1.*delta,2)*pow(mrho,8) +
             pow(mpion,2)*pow(mrho,4)*(8.000000000000002*pow(2. - 1.*delta,2)*pow(mrho,2) +
                (-32. + 31.999999999999996*delta - 8.*pow(delta,2))*s))*log(abs(-2.*pow(mpion,2) + s + tmin)))/
         (pow(mrho,4)*(-1.*pow(mpion,2) + 1.*s)) +
        (0.25*(0. + 8.*pow(mpion,2)*pow(mrho,6) + 4.*pow(mrho,8) - 8.*pow(mrho,6)*s +
             delta*pow(mrho,4)*(8.*pow(mpion,4) - 8.*pow(mrho,4) + pow(mpion,2)*(8.*pow(mrho,2) - 16.*s) +
                8.*pow(mrho,2)*s + 8.*pow(s,2)) +
             pow(delta,2)*pow(mrho,4)*(-4.*pow(mpion,4) + 3.*pow(mrho,4) - 2.*pow(mrho,2)*s -
                4.*pow(s,2) + pow(mpion,2)*(-6.*pow(mrho,2) + 8.*s)))*log(abs(-2.*pow(mpion,2) + s + tmin)))/
         pow(mrho,6) - (0.5*(0. + pow(mpion,2)*(4.*pow(mrho,6) - 8.*C4*pow(mrho,8)) - 4.*pow(mrho,6)*s +
             8.*C4*pow(mrho,8)*s + pow(delta,2)*pow(mrho,4)*
              (-2.*pow(mpion,4) + 1.*pow(mrho,4) - 2.*pow(s,2) + pow(mpion,2)*(-4.*pow(mrho,2) + 4.*s)) +
             delta*pow(mrho,4)*(4.*pow(mpion,4) +
                pow(mpion,2)*(6.*pow(mrho,2) + 4.*C4*pow(mrho,4) - 8.*s) + 2.*pow(mrho,2)*s + 4.*pow(s,2) +
                pow(mrho,4)*(-2. - 4.*C4*s)))*log(abs(-2.*pow(mpion,2) + s + tmin)))/pow(mrho,6)))/
    (16.*Pi*(pow(mpion,4) + pow(pow(mrho,2) - s,2) - 2*pow(mpion,2)*(pow(mrho,2) + s)));

  return sigma / spin_deg_factor;
}

//C14
double total_xsection_C14(double s) {
  double sigma;
  double spin_deg_factor = 3.0;
  double tmin = min_mandelstam_t(s, mpion, mrho, mpion);
  double tmax = max_mandelstam_t(s, mpion, mrho, mpion);

  sigma = (pow(Const,2)*pow(g_POR,4)*((-0.125*pow(momega,8) - 0.125*pow(mpion,8) + 0.25*pow(mpion,6)*pow(mrho,2) - 0.125*pow(mpion,4)*pow(mrho,4) +
	          pow(momega,6)*(0.5*pow(mpion,2) + 0.25*pow(mrho,2) - 0.25*s) + pow(momega,2)*pow(mpion,2)*(0.25*pow(mrho,4) + 0.25*pow(mpion,2)*s - 0.25*pow(mrho,2)*s) +
	          pow(momega,4)*(-0.5*pow(mpion,4) - 0.125*pow(mrho,4) + pow(mpion,2)*(-0.5*pow(mrho,2) + 0.5*s) + 0.25*pow(mrho,2)*s - 0.25*pow(s,2)))/(1.*pow(momega,2) - 1.*tmin)
	         - 0.125*(3.*pow(momega,4) + 4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2) +
	          pow(momega,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s))*tmin - 0.125*(pow(momega,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s)*pow(tmin,2) -
	       0.041666666666666664*pow(tmin,3) + (0.125*pow(momega,8) + 0.125*pow(mpion,8) - 0.25*pow(mpion,6)*pow(mrho,2) + 0.125*pow(mpion,4)*pow(mrho,4) +
	          pow(momega,6)*(-0.5*pow(mpion,2) - 0.25*pow(mrho,2) + 0.25*s) + pow(momega,2)*pow(mpion,2)*(-0.25*pow(mrho,4) - 0.25*pow(mpion,2)*s + 0.25*pow(mrho,2)*s) +
	          pow(momega,4)*(0.5*pow(mpion,4) + 0.125*pow(mrho,4) + pow(mpion,2)*(0.5*pow(mrho,2) - 0.5*s) - 0.25*pow(mrho,2)*s + 0.25*pow(s,2)))/(1.*pow(momega,2) - 1.*tmax)\
	        + 0.125*(3.*pow(momega,4) + 4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2) +
	          pow(momega,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s))*tmax + 0.125*(pow(momega,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s)*pow(tmax,2) +
	       0.041666666666666664*pow(tmax,3) - 0.25*(2.*pow(momega,6) + pow(momega,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
	          pow(mpion,2)*(-1.*pow(mrho,4) - 1.*pow(mpion,2)*s + pow(mrho,2)*s) +
	          pow(momega,2)*(4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2)))*log(fabs(-1.*pow(momega,2) + tmin)) -
	       (HeavisideTheta(-momega + sqrt(s))*(tmin*(0.125*pow(mpion,8) - 0.25*pow(mpion,6)*pow(mrho,2) +
	               pow(mpion,4)*(0.125*pow(mrho,4) + pow(momega,2)*(0.25*pow(mrho,2) - 0.5*s) - 0.25*pow(mrho,2)*s + s*(1.*s - 0.125*tmin)) +
	               pow(mpion,2)*s*(1.*pow(momega,4) - 0.25*pow(mrho,4) + s*(-1.5*s - 0.75*tmin) + pow(mrho,2)*(1.5*s + 0.125*tmin) +
	                  pow(momega,2)*(-1.0000000000000002*pow(mrho,2) + 0.5*tmin)) +
	               s*(-0.25*pow(momega,6) + pow(momega,4)*(0.25*pow(mrho,2) - 0.5*s - 0.125*tmin) +
	                  pow(momega,2)*(0.5*pow(s,2) - 0.25*s*tmin + (0.125*pow(mrho,2) - 0.08333333333333333*tmin)*tmin) +
	                  s*(0.125*pow(mrho,4) + 0.375*pow(s,2) + pow(mrho,2)*(-0.5*s - 0.25*tmin) + 0.5*s*tmin + 0.16666666666666666*pow(tmin,2)))) +
	            (-0.25*pow(momega,8)*s + pow(momega,6)*(1.*pow(mpion,2) + 0.25*pow(mrho,2) - 0.5*s)*s + pow(mpion,4)*s*(0.25*pow(mpion,4) - 0.25*pow(mrho,2)*s) +
	               pow(momega,4)*(pow(mpion,4)*(0.25*pow(mrho,2) - 0.5*s) - 1.0000000000000002*pow(mpion,2)*pow(mrho,2)*s + 0.5*pow(s,3)) +
	               pow(momega,2)*(-0.25*pow(mpion,8) + 0.5*pow(mpion,4)*pow(s,2) + pow(mpion,2)*(1.*pow(mrho,2) - 1.*s)*pow(s,2) + (-0.25*pow(mrho,2) + 0.25*s)*pow(s,3)))*
	             log(fabs(-1.*pow(momega,2) + tmin))))/pow(pow(momega,2) - 1.*s,2) +
	       0.25*(2.*pow(momega,6) + pow(momega,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) + pow(mpion,2)*(-1.*pow(mrho,4) - 1.*pow(mpion,2)*s + pow(mrho,2)*s) +
	          pow(momega,2)*(4.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2)))*log(fabs(-1.*pow(momega,2) + tmax)) +
	       (HeavisideTheta(-momega + sqrt(s))*(tmax*(0.125*pow(mpion,8) - 0.25*pow(mpion,6)*pow(mrho,2) +
	               pow(mpion,4)*(0.125*pow(mrho,4) + pow(momega,2)*(0.25*pow(mrho,2) - 0.5*s) - 0.25*pow(mrho,2)*s + s*(1.*s - 0.125*tmax)) +
	               pow(mpion,2)*s*(1.*pow(momega,4) - 0.25*pow(mrho,4) + s*(-1.5*s - 0.75*tmax) + pow(mrho,2)*(1.5*s + 0.125*tmax) +
	                  pow(momega,2)*(-1.0000000000000002*pow(mrho,2) + 0.5*tmax)) +
	               s*(-0.25*pow(momega,6) + pow(momega,4)*(0.25*pow(mrho,2) - 0.5*s - 0.125*tmax) +
	                  pow(momega,2)*(0.5*pow(s,2) - 0.25*s*tmax + (0.125*pow(mrho,2) - 0.08333333333333333*tmax)*tmax) +
	                  s*(0.125*pow(mrho,4) + 0.375*pow(s,2) + pow(mrho,2)*(-0.5*s - 0.25*tmax) + 0.5*s*tmax + 0.16666666666666666*pow(tmax,2)))) +
	            (-0.25*pow(momega,8)*s + pow(momega,6)*(1.*pow(mpion,2) + 0.25*pow(mrho,2) - 0.5*s)*s + pow(mpion,4)*s*(0.25*pow(mpion,4) - 0.25*pow(mrho,2)*s) +
	               pow(momega,4)*(pow(mpion,4)*(0.25*pow(mrho,2) - 0.5*s) - 1.0000000000000002*pow(mpion,2)*pow(mrho,2)*s + 0.5*pow(s,3)) +
	               pow(momega,2)*(-0.25*pow(mpion,8) + 0.5*pow(mpion,4)*pow(s,2) + pow(mpion,2)*(1.*pow(mrho,2) - 1.*s)*pow(s,2) + (-0.25*pow(mrho,2) + 0.25*s)*pow(s,3)))*
	             log(fabs(-1.*pow(momega,2) + tmax))))/pow(pow(momega,2) - 1.*s,2)))/(16.*Pi*(pow(mpion,4) + pow(pow(mrho,2) - s,2) - 2*pow(mpion,2)*(pow(mrho,2) + s)));

  return sigma / spin_deg_factor;
}

//C15
double total_xsection_C15(double s) {
  double sigma;
  double spin_deg_factor = 3.0;
  double tmin = min_mandelstam_t(s, mpion, mrho, mpion);
  double tmax = max_mandelstam_t(s, mpion, mrho, mpion);

  sigma = HeavisideTheta(-momega + sqrt(s))*((0.0024867959858108648*pow(Const,2)*pow(g_POR,4)*(pow(mpion,8)*(1.*tmax - 1.*tmin) + pow(mpion,6)*pow(mrho,2)*(-2.*tmax + 2.*tmin) +
       pow(mpion,4)*(pow(mrho,4)*(1.*tmax - 1.*tmin) + s*(4.*s*tmax - 1.*pow(tmax,2) - 4.*s*tmin + 1.*pow(tmin,2))) +
       pow(s,2)*(1.*pow(s,2)*tmax + 1.*s*pow(tmax,2) + 0.6666666666666666*pow(tmax,3) + pow(mrho,4)*(1.*tmax - 1.*tmin) -
          1.*pow(s,2)*tmin - 1.*s*pow(tmin,2) - 0.6666666666666666*pow(tmin,3) +
          pow(mrho,2)*(-2.*s*tmax - 1.*pow(tmax,2) + 2.*s*tmin + 1.*pow(tmin,2))) +
       pow(mpion,2)*s*(pow(mrho,4)*(-2.*tmax + 2.*tmin) + pow(mrho,2)*(4.*s*tmax + 1.*pow(tmax,2) - 4.*s*tmin - 1.*pow(tmin,2)) +
          s*(-4.*s*tmax - 2.*pow(tmax,2) + 4.*s*tmin + 2.*pow(tmin,2)))))/
   ((pow(pow(momega,2) - 1.*s,2)*(pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-2.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s + pow(s,2)))));

  return sigma / spin_deg_factor;
}

//C16
double total_xsection_C16(double s) {
  double sigma;
  double spin_deg_factor = 3.0;
  double tmin = min_mandelstam_t(s, mpion, mrho, mpion);
  double tmax = max_mandelstam_t(s, mpion, mrho, mpion);

  sigma = (0.0024868*pow(Const,2)*pow(g_POR,4)*((pow(momega,8) + pow(mpion,4)*pow(pow(mpion,2) - pow(mrho,2),2) -
          2*pow(momega,6)*(2*pow(mpion,2) + pow(mrho,2) - s) -
          2*pow(momega,2)*pow(mpion,2)*(pow(mrho,4) + pow(mpion,2)*s - pow(mrho,2)*s) +
          pow(momega,4)*(4*pow(mpion,4) + pow(mrho,4) + 4*pow(mpion,2)*(pow(mrho,2) - s) - 2*pow(mrho,2)*s + 2*pow(s,2)))/
        (pow(momega,2) - tmax) + 3*pow(momega,4)*tmax - 8*pow(momega,2)*pow(mpion,2)*tmax + 4*pow(mpion,4)*tmax -
       4*pow(momega,2)*pow(mrho,2)*tmax + 4*pow(mpion,2)*pow(mrho,2)*tmax + pow(mrho,4)*tmax + 4*pow(momega,2)*s*tmax -
       4*pow(mpion,2)*s*tmax - 2*pow(mrho,2)*s*tmax + 2*pow(s,2)*tmax + pow(momega,2)*pow(tmax,2) - 2*pow(mpion,2)*pow(tmax,2) -
       pow(mrho,2)*pow(tmax,2) + s*pow(tmax,2) + pow(tmax,3)/3. -
       (pow(momega,8) + pow(mpion,4)*pow(pow(mpion,2) - pow(mrho,2),2) - 2*pow(momega,6)*(2*pow(mpion,2) + pow(mrho,2) - s) -
          2*pow(momega,2)*pow(mpion,2)*(pow(mrho,4) + pow(mpion,2)*s - pow(mrho,2)*s) +
          pow(momega,4)*(4*pow(mpion,4) + pow(mrho,4) + 4*pow(mpion,2)*(pow(mrho,2) - s) - 2*pow(mrho,2)*s + 2*pow(s,2)))/
        (pow(momega,2) - tmin) - 3*pow(momega,4)*tmin + 8*pow(momega,2)*pow(mpion,2)*tmin - 4*pow(mpion,4)*tmin +
       4*pow(momega,2)*pow(mrho,2)*tmin - 4*pow(mpion,2)*pow(mrho,2)*tmin - pow(mrho,4)*tmin - 4*pow(momega,2)*s*tmin +
       4*pow(mpion,2)*s*tmin + 2*pow(mrho,2)*s*tmin - 2*pow(s,2)*tmin - pow(momega,2)*pow(tmin,2) + 2*pow(mpion,2)*pow(tmin,2) +
       pow(mrho,2)*pow(tmin,2) - s*pow(tmin,2) - pow(tmin,3)/3. +
       2*(2*pow(momega,6) - 3*pow(momega,4)*(2*pow(mpion,2) + pow(mrho,2) - s) -
          pow(mpion,2)*(pow(mrho,4) + pow(mpion,2)*s - pow(mrho,2)*s) +
          pow(momega,2)*(4*pow(mpion,4) + pow(mrho,4) + 4*pow(mpion,2)*(pow(mrho,2) - s) - 2*pow(mrho,2)*s + 2*pow(s,2)))*
        log(fabs(-pow(momega,2) + tmax)) - 2*(2*pow(momega,6) - 3*pow(momega,4)*(2*pow(mpion,2) + pow(mrho,2) - s) -
          pow(mpion,2)*(pow(mrho,4) + pow(mpion,2)*s - pow(mrho,2)*s) +
          pow(momega,2)*(4*pow(mpion,4) + pow(mrho,4) + 4*pow(mpion,2)*(pow(mrho,2) - s) - 2*pow(mrho,2)*s + 2*pow(s,2)))*
        log(fabs(-pow(momega,2) + tmin))))/((pow(mpion,4) + pow(pow(mrho,2) - s,2) - 2*pow(mpion,2)*(pow(mrho,2) + s)));

  return sigma / spin_deg_factor;
}

//C21
double total_xsection_C21(double s) {
  double sigma;
  double spin_deg_factor = 1.0;
  double tmin = min_mandelstam_t(s, mpion, mpion, mrho);
  double tmax = max_mandelstam_t(s, mpion, mpion, mrho);

  sigma = (-(pow(Const,2)*pow(ghat,4)*(0. + (0.03125*pow(eta1 -
							1.*eta2,2)*(eta1*eta2*
							              (-2.*pow(ma1,8) - 2.*pow(mpion,8) +
							2.*pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) +
							                pow(ma1,2)*pow(mpion,2)*(8.*pow(mpion,4) -
							8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
							                pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) +
							8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
							             pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) -
							2.*pow(mpion,6)*pow(mrho,2) + 1.*pow(mpion,4)*pow(mrho,4) +
							                pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) +
							2.*s) +
							                pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) +
							pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
							2.*pow(s,2)) +
							                pow(ma1,2)*(-4.*pow(mpion,6) -
							2.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(6.*pow(mrho,2) +
							2.*s))) +
							             pow(eta1,2)*(1.*pow(ma1,8) +
							pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
							                pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) +
							pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
							2.*pow(s,2)) +
							                pow(ma1,2)*(-4.*pow(mpion,6) +
							2.*pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(2.*pow(mrho,2) -
							2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
							                pow(mpion,2)*(1.*pow(mpion,6) +
							2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s
							+ pow(mpion,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s)))))/
							         (1.*pow(ma1,2) - 1.*tmin) + (1.*pow(-2. +
							delta,2)*pow(mpion,2)*(1.*pow(mpion,2) -
							0.25*pow(mrho,2)))/(1.*pow(mpion,2) - 1.*tmin) +
							        (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) -
							0.25*pow(mrho,2)))/(1.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s -
							1.*tmin) -
							        (0.5*pow(-2. + delta,2)*pow(mpion,2)*tmin)/pow(mrho,2) -
							0.25*(-2. + delta)*(eta1 - 1.*eta2)*
							         (-0.5*eta2*pow(ma1,2) + 1.*eta1*pow(mpion,2) +
							0.5*eta2*pow(mrho,2) + 0.5*eta1*s - 1.*eta2*s)*tmin +
							        (0.25*(pow(mpion,2)*(12. + 1.*pow(delta,2) -
							16.*C4*pow(mrho,2) + delta*(-8. + 8.*C4*pow(mrho,2))) +
							             (-4. - 3.*pow(delta,2) - 16.*C4*pow(mrho,2) +
							delta*(8. + 8.*C4*pow(mrho,2)))*s)*tmin)/pow(mrho,2) -
							        0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(pow(ma1,2) -
							1.*s) + eta1*(-2.*pow(mpion,2) + s))*tmin -
							        0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-1.*pow(ma1,2) -
							2.*pow(mpion,2) - 1.*s) + eta1*(2.*pow(mpion,2) + s))*tmin +
							        0.25*(-2. + delta)*(eta1 - 1.*eta2)*(eta1*(1.*pow(mpion,2) -
							0.5*s) + eta2*(-0.5*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2)
							+ 1.*s))*tmin +
							        (0.25*(-2. + 1.*delta)*(-8.*C4*pow(mrho,4) +
							pow(mpion,2)*(2. + 1.*delta - 8.*C4*pow(mrho,2)) + (-2. -
							3.*delta)*s + pow(mrho,2)*(2. + 1.*delta + 16.*C4*s))*tmin)/
							         pow(mrho,2) + (0.25*(32*pow(C4,2)*pow(mrho,8) +
							2*pow(delta,2)*pow(s,2) + 8*C4*pow(mrho,6)*(-6 + delta - 8*C4*s) +
							2*delta*pow(mrho,2)*s*(-6 + delta - 8*C4*s) +
							             pow(mrho,4)*(12 - pow(delta,2) + 8*C4*(6 + delta)*s +
							32*pow(C4,2)*pow(s,2)))*tmin)/pow(mrho,4) -
							        (1.*(1.*eta1 - 1.*eta2)*(eta2*(0.75*pow(mrho,4) -
							0.125*delta*pow(mrho,4) - 1.*C4*pow(mrho,6) +
							pow(ma1,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) +
							                pow(mpion,2)*(-1.*pow(mrho,2) +
							2.*C4*pow(mrho,4)) - 0.25*pow(mrho,2)*s -
							0.375*delta*pow(mrho,2)*s + 2.*C4*pow(mrho,4)*s +
							0.25*delta*pow(s,2) -
							                1.*C4*pow(mrho,2)*pow(s,2)) +
							eta1*(0.5*pow(mrho,4) - 1.*C4*pow(mrho,6) +
							pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
							                pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4))
							- 0.25*delta*pow(s,2) +
							1.*C4*pow(mrho,2)*pow(s,2)))*tmin)/pow(mrho,2) +
							        0.0625*pow(eta1 -
							1.*eta2,2)*(pow(eta2,2)*(pow(Gammaa1,2)*pow(ma1,2) -
							1.*pow(ma1,4) - 2.*pow(mpion,4) - 2.*pow(mpion,2)*pow(mrho,2) +
							2.*pow(mrho,4) +
							              pow(ma1,2)*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s) -
							3.*pow(mrho,2)*s + pow(s,2)) +
							           pow(eta1,2)*(pow(Gammaa1,2)*pow(ma1,2) -
							1.*pow(ma1,4) - 2.*pow(mpion,4) - 2.*pow(mpion,2)*pow(mrho,2) +
							              pow(ma1,2)*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s) +
							pow(mrho,2)*s + pow(s,2)) +
							           eta1*eta2*(-2.*pow(Gammaa1,2)*pow(ma1,2) +
							2.*pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) +
							2.*pow(mrho,4) - 2.*pow(s,2) +
							              pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) +
							2.*s)))*tmin +
							        0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-6.*pow(ma1,4) -
							12.*pow(mpion,4) + 2.*pow(mrho,4) + pow(ma1,2)*(16.*pow(mpion,2)
							- 8.*s) + 8.*pow(mpion,2)*s +
							              4.*pow(mrho,2)*s - 4.*pow(s,2)) +
							pow(eta1,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) +
							pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
							              2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) -
							4.*pow(mrho,2) + 4.*s)) +
							           pow(eta2,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) +
							pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) -
							2.*pow(mrho,2)*s + 2.*pow(s,2) +
							              pow(ma1,2)*(-8.*pow(mpion,2) + 4.*pow(mrho,2) +
							4.*s)))*tmin - (0.125*(-2. + 1.*delta)*(2. + 1.*delta -
							8.*C4*pow(mrho,2))*pow(tmin,2))/pow(mrho,2) -
							        0.5*pow(1.*eta1 - 1.*eta2,2)*(-0.5 +
							1.*C4*pow(mrho,2))*pow(tmin,2) -
							        (1.*(0.5 - 0.125*pow(delta,2) - 2.*C4*pow(mrho,2) +
							1.*C4*delta*pow(mrho,2))*pow(tmin,2))/pow(mrho,2) +
							        0.0625*pow(1.*eta1 - 1.*eta2,4)*(1.*pow(mpion,2) +
							0.5*pow(mrho,2) - 0.5*s)*pow(tmin,2) +
							        0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) +
							2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) + eta1*(pow(ma1,2) -
							2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*
							         pow(tmin,2) + 0.010416666666666666*pow(eta1 -
							1.*eta2,4)*pow(tmin,3) - 0.020833333333333332*pow(1.*eta1 -
							1.*eta2,4)*pow(tmin,3) +
							        0.03125*pow(eta1 -
							1.*eta2,2)*(eta1*eta2*(2.*pow(Gammaa1,2)*pow(ma1,2) -
							6.*pow(ma1,4) - 4.*pow(mpion,4) + 2.*pow(mrho,4) +
							              pow(ma1,2)*(8.*pow(mpion,2) - 8.*s) +
							4.*pow(mrho,2)*s - 4.*pow(s,2)) +
							           pow(eta1,2)*(-1.*pow(Gammaa1,2)*pow(ma1,2) +
							3.*pow(ma1,4) + 2.*pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) +
							pow(mrho,4) - 4.*pow(mrho,2)*s +
							              2.*pow(s,2) + pow(ma1,2)*(-4.*pow(mpion,2) -
							4.*pow(mrho,2) + 4.*s)) +
							           pow(eta2,2)*(-1.*pow(Gammaa1,2)*pow(ma1,2) +
							3.*pow(ma1,4) + 2.*pow(mpion,4) - 2.*pow(mpion,2)*pow(mrho,2) +
							pow(mrho,4) - 2.*pow(mrho,2)*s +
							              2.*pow(s,2) + pow(ma1,2)*(-4.*pow(mpion,2) +
							4.*pow(mrho,2) + 4.*s)))*(-1.*pow(mrho,2) + s + tmin) -
							        0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) -
							1.*pow(mrho,2) - 1.*s) + eta1*(pow(ma1,2) - 1.*pow(mrho,2) +
							s))*pow(-1.*pow(mrho,2) + s + tmin,2) +
							        0.010416666666666666*pow(eta1 -
							1.*eta2,4)*pow(-1.*pow(mrho,2) + s + tmin,3) +
							        0.25*(eta1 - 1.*eta2)*(1.*eta1 - 1.*eta2)*(-1. +
							2.*C4*pow(mrho,2))*pow(pow(ma1,2) - 2.*pow(mpion,2) -
							1.*pow(mrho,2) + s + tmin,2) -
							        (2.*(1.*eta1 - 1.*eta2)*(eta2*(0.375*pow(mrho,4) -
							0.0625*delta*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
							pow(ma1,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) +
							                pow(mpion,2)*(-0.5*pow(mrho,2) +
							1.*C4*pow(mrho,4)) - 0.125*pow(mrho,2)*s -
							0.1875*delta*pow(mrho,2)*s + 1.*C4*pow(mrho,4)*s +
							0.125*delta*pow(s,2) -
							                0.5*C4*pow(mrho,2)*pow(s,2)) +
							eta1*(0.25*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
							pow(mpion,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) +
							                pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4))
							- 0.125*delta*pow(s,2) + 0.5*C4*pow(mrho,2)*pow(s,2)))*
							           (1.*pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) +
							1.*s + 1.*tmin))/pow(mrho,2) +
							        (2.*(1.*eta1 - 1.*eta2)*Gammaa1*ma1*(eta2*(0.375*pow(mrho,4)
							- 0.0625*delta*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
							pow(ma1,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) +
							                pow(mpion,2)*(-0.5*pow(mrho,2) +
							1.*C4*pow(mrho,4)) - 0.125*pow(mrho,2)*s -
							0.1875*delta*pow(mrho,2)*s + 1.*C4*pow(mrho,4)*s +
							0.125*delta*pow(s,2) -
							                0.5*C4*pow(mrho,2)*pow(s,2)) +
							eta1*(0.25*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
							pow(mpion,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) +
							                pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4))
							- 0.125*delta*pow(s,2) + 0.5*C4*pow(mrho,2)*pow(s,2)))*
							           atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2)
							+ s + tmin)/(Gammaa1*ma1)))/pow(mrho,2) +
							        (0.25*(-2. + delta)*(eta1 -
							1.*eta2)*Gammaa1*ma1*(eta2*(-1.*pow(ma1,6) +
							pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(ma1,2) + 0.5*pow(mrho,2) -
							1.*s) +
							                pow(ma1,4)*(2.*pow(mpion,2) + 0.5*pow(mrho,2) -
							1.*s) + pow(mpion,4)*(-1.5*pow(mrho,2) + 1.*s) +
							                pow(ma1,2)*pow(mpion,2)*(-1.*pow(mpion,2) -
							1.*pow(mrho,2) + 2.*s)) +
							             eta1*(pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,2) +
							0.5*s) + pow(ma1,4)*(1.*pow(mpion,2) + 0.5*s) +
							pow(ma1,2)*(-2.*pow(mpion,4) - 1.*pow(mpion,2)*s) +
							                pow(mpion,2)*(1.*pow(mpion,4) - 1.*pow(mrho,4) +
							pow(mpion,2)*(2.*pow(mrho,2) - 1.5*s) + 1.*pow(mrho,2)*s)))*
							           atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2)
							+ s + tmin)/(Gammaa1*ma1)))/
							         (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) -
							2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) -
							        (0.25*(-2. + delta)*(eta1 -
							1.*eta2)*Gammaa1*ma1*(eta2*(-1.*pow(ma1,6) -
							2.*pow(mpion,4)*pow(mrho,2) - 1.*pow(mpion,2)*pow(mrho,4) +
							                pow(ma1,4)*(2.*pow(mpion,2) + 2.*pow(mrho,2) -
							1.5*s) + 2.5*pow(mpion,4)*s + 3.*pow(mpion,2)*pow(mrho,2)*s +
							0.5*pow(mrho,4)*s -
							                2.*pow(mpion,2)*pow(s,2) -
							1.*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) +
							pow(Gammaa1,2)*(-1.*pow(ma1,4) + 0.5*pow(ma1,2)*s) +
							                pow(ma1,2)*(-1.*pow(mpion,4) - 1.*pow(mrho,4) +
							1.*pow(mrho,2)*s + pow(mpion,2)*(-2.*pow(mrho,2) + 1.*s))) +
							             eta1*(1.*pow(mpion,6) + 4.*pow(mpion,4)*pow(mrho,2)
							+ 1.*pow(mpion,2)*pow(mrho,4) +
							pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,2) - 0.5*s) +
							                pow(ma1,4)*(1.*pow(mpion,2) - 0.5*s) -
							4.5*pow(mpion,4)*s - 4.*pow(mpion,2)*pow(mrho,2)*s -
							0.5*pow(mrho,4)*s + 3.*pow(mpion,2)*pow(s,2) +
							                1.*pow(mrho,2)*pow(s,2) - 0.5*pow(s,3) +
							pow(ma1,2)*(-2.*pow(mpion,4) + (1.*pow(mrho,2) - 1.*s)*s +
							pow(mpion,2)*(-2.*pow(mrho,2) + 3.*s))))*
							           atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2)
							+ s + tmin)/(Gammaa1*ma1)))/
							         (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4)
							+ 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 2.*pow(mpion,2)*s
							- 2.*pow(mrho,2)*s + pow(s,2) +
							           pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) +
							        (0.03125*pow(eta1 -
							1.*eta2,2)*(pow(eta2,2)*(pow(Gammaa1,4)*pow(ma1,4) + pow(ma1,8)
							+ pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) +
							                pow(mpion,4)*pow(mrho,4) +
							pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
							                pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) +
							pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
							2.*pow(s,2)) +
							                pow(ma1,2)*(-4.*pow(mpion,6) -
							2.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(6.*pow(mrho,2) +
							2.*s)) +
							                pow(Gammaa1,2)*pow(ma1,2)*(-6.*pow(ma1,4) -
							6.*pow(mpion,4) - 1.*pow(mrho,4) + pow(ma1,2)*(12.*pow(mpion,2)
							- 6.*pow(mrho,2) - 6.*s) +
							                   2.*pow(mrho,2)*s - 2.*pow(s,2) +
							pow(mpion,2)*(6.*pow(mrho,2) + 4.*s))) +
							             eta1*eta2*(-2.*pow(Gammaa1,4)*pow(ma1,4) -
							2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) +
							pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) +
							                pow(ma1,2)*pow(mpion,2)*(8.*pow(mpion,4) -
							8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
							                pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) +
							8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2)) +
							                pow(Gammaa1,2)*pow(ma1,2)*(12.*pow(ma1,4) +
							12.*pow(mpion,4) - 2.*pow(mrho,4) - 8.*pow(mpion,2)*s -
							4.*pow(mrho,2)*s + 4.*pow(s,2) +
							                   pow(ma1,2)*(-24.*pow(mpion,2) + 12.*s))) +
							pow(eta1,2)*
							              (pow(Gammaa1,4)*pow(ma1,4) + pow(ma1,8) +
							pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
							                pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) +
							pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
							2.*pow(s,2)) +
							                pow(ma1,2)*(-4.*pow(mpion,6) +
							2.*pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(2.*pow(mrho,2) -
							2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
							                pow(Gammaa1,2)*pow(ma1,2)*(-6.*pow(ma1,4) -
							6.*pow(mpion,4) - 1.*pow(mrho,4) + pow(ma1,2)*(12.*pow(mpion,2)
							+ 6.*pow(mrho,2) - 6.*s) +
							                   4.*pow(mrho,2)*s - 2.*pow(s,2) +
							pow(mpion,2)*(-6.*pow(mrho,2) + 4.*s)) +
							                pow(mpion,2)*(pow(mpion,6) +
							2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s
							+ pow(mpion,2)*(pow(mrho,4) - 2.*pow(mrho,2)*s))))*
							           atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2)
							+ s + tmin)/(Gammaa1*ma1)))/(Gammaa1*ma1) -
							        (0.0625*pow(eta1 -
							1.*eta2,2)*Gammaa1*ma1*(eta1*eta2*(-2.*pow(Gammaa1,4)*pow(ma1,4) +
							14.*pow(ma1,8) + 14.*pow(mpion,8) + 28.*pow(mpion,6)*pow(mrho,2) +
							                20.*pow(mpion,4)*pow(mrho,4) +
							10.*pow(mpion,2)*pow(mrho,6) + 2.*pow(mrho,8) -
							16.*pow(mpion,6)*s - 16.*pow(mpion,4)*pow(mrho,2)*s -
							                12.*pow(mpion,2)*pow(mrho,4)*s - 4.*pow(mrho,6)*s
							- 4.*pow(mpion,4)*pow(s,2) -
							6.*pow(mpion,2)*pow(mrho,2)*pow(s,2) + 8.*pow(mpion,2)*pow(s,3) +
							                4.*pow(mrho,2)*pow(s,3) - 2.*pow(s,4) +
							pow(ma1,6)*(-56.*pow(mpion,2) - 28.*pow(mrho,2) + 28.*s) +
							                pow(ma1,4)*(84.*pow(mpion,4) + 24.*pow(mrho,4) +
							pow(mpion,2)*(84.*pow(mrho,2) - 72.*s) - 36.*pow(mrho,2)*s +
							12.*pow(s,2)) +
							                pow(Gammaa1,2)*pow(ma1,2)*(-4.*pow(ma1,4) -
							4.*pow(mpion,4) + pow(ma1,2)*(8.*pow(mpion,2) + 4.*pow(mrho,2) -
							4.*s) + (4.*pow(mrho,2) - 4.*s)*s +
							                   pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s)) +
							pow(ma1,2)*
							                 (-56.*pow(mpion,6) - 10.*pow(mrho,6) +
							18.*pow(mrho,4)*s - 6.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3) +
							pow(mpion,4)*(-84.*pow(mrho,2) + 60.*s) +
							                   pow(mpion,2)*(-48.*pow(mrho,4) +
							60.*pow(mrho,2)*s - 12.*pow(s,2)))) +
							             pow(eta1,2)*(1.*pow(Gammaa1,4)*pow(ma1,4) -
							7.*pow(ma1,8) - 7.*pow(mpion,8) - 14.*pow(mpion,6)*pow(mrho,2) -
							7.*pow(mpion,4)*pow(mrho,4) -
							                2.*pow(mpion,2)*pow(mrho,6) +
							pow(ma1,6)*(28.*pow(mpion,2) + 14.*pow(mrho,2) - 14.*s) +
							8.*pow(mpion,6)*s + 11.*pow(mpion,4)*pow(mrho,2)*s +
							                6.*pow(mpion,2)*pow(mrho,4)*s + 1.*pow(mrho,6)*s
							+ 2.*pow(mpion,4)*pow(s,2) - 1.*pow(mrho,4)*pow(s,2) -
							4.*pow(mpion,2)*pow(s,3) -
							                1.*pow(mrho,2)*pow(s,3) + 1.*pow(s,4) +
							pow(Gammaa1,2)*pow(ma1,2)*
							                 (2.*pow(ma1,4) + 2.*pow(mpion,4) +
							1.*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 4.*s) -
							1.*pow(mrho,2)*s + 2.*pow(s,2) +
							                   pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2)
							+ 2.*s)) +
							                pow(ma1,4)*(-42.*pow(mpion,4) - 9.*pow(mrho,4) +
							21.*pow(mrho,2)*s - 6.*pow(s,2) + pow(mpion,2)*(-42.*pow(mrho,2)
							+ 36.*s)) +
							                pow(ma1,2)*(28.*pow(mpion,6) + 2.*pow(mrho,6) +
							pow(mpion,4)*(42.*pow(mrho,2) - 30.*s) - 9.*pow(mrho,4)*s +
							6.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) +
							                   pow(mpion,2)*(18.*pow(mrho,4) -
							36.*pow(mrho,2)*s + 6.*pow(s,2)))) +
							             pow(eta2,2)*(1.*pow(Gammaa1,4)*pow(ma1,4) -
							7.*pow(ma1,8) - 7.*pow(mpion,8) - 14.*pow(mpion,6)*pow(mrho,2) -
							1.*pow(mpion,4)*pow(mrho,4) +
							                6.*pow(mpion,2)*pow(mrho,6) + 2.*pow(mrho,8) +
							pow(ma1,6)*(28.*pow(mpion,2) + 14.*pow(mrho,2) - 14.*s) +
							8.*pow(mpion,6)*s -
							                1.*pow(mpion,4)*pow(mrho,2)*s -
							16.*pow(mpion,2)*pow(mrho,4)*s - 7.*pow(mrho,6)*s +
							2.*pow(mpion,4)*pow(s,2) +
							                14.*pow(mpion,2)*pow(mrho,2)*pow(s,2) +
							9.*pow(mrho,4)*pow(s,2) - 4.*pow(mpion,2)*pow(s,3) -
							5.*pow(mrho,2)*pow(s,3) + 1.*pow(s,4) +
							                pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,4) +
							2.*pow(mpion,4) + 3.*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2)
							- 4.*s) - 5.*pow(mrho,2)*s + 2.*pow(s,2) +
							                   pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2)
							+ 2.*s)) +
							                pow(ma1,4)*(-42.*pow(mpion,4) - 3.*pow(mrho,4) +
							9.*pow(mrho,2)*s - 6.*pow(s,2) + pow(mpion,2)*(-42.*pow(mrho,2)
							+ 36.*s)) +
							                pow(ma1,2)*(28.*pow(mpion,6) - 4.*pow(mrho,6) +
							pow(mpion,4)*(42.*pow(mrho,2) - 30.*s) + 9.*pow(mrho,4)*s -
							6.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) +
							                   pow(mpion,2)*(6.*pow(mrho,4) -
							12.*pow(mrho,2)*s + 6.*pow(s,2)))))*atan((pow(ma1,2) -
							2.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmin)/(Gammaa1*ma1)))/
							         (pow(Gammaa1,2)*pow(ma1,2) + 4.*pow(ma1,4) +
							4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
							4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s +
							           pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) -
							4.*pow(mrho,2) + 4.*s)) +
							        0.0625*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-4.*pow(ma1,6) +
							pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) +
							              pow(mpion,2)*(4.*pow(mpion,4) - 4.*pow(mrho,4) -
							2.*pow(mpion,2)*s + 2.*pow(mrho,2)*s) +
							              pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) +
							8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
							           pow(eta1,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) +
							pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(pow(mrho,2) - 1.*s)*s
							+ pow(mpion,4)*(-3.*pow(mrho,2) + s) +
							              pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) +
							3.*s) +
							              pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) +
							pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
							2.*pow(s,2))) +
							           pow(eta2,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) -
							1.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(3.*pow(mrho,2) + s) +
							              pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) +
							3.*s) +
							              pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) +
							pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
							2.*pow(s,2))))*log(fabs(-1.*pow(ma1,2) + tmin)) -
							        (0.25*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-0.5*pow(ma1,6) -
							0.5*pow(mpion,6) + 0.5*pow(mpion,4)*pow(mrho,2) +
							                pow(ma1,4)*(0.5*pow(mpion,2) + 0.5*pow(mrho,2) -
							1.*s) + pow(ma1,2)*pow(mpion,2)*(0.5*pow(mpion,2) +
							1.*pow(mrho,2) - 1.*s)) +
							             eta1*(pow(ma1,4)*(1.*pow(mpion,2) + 0.5*s) +
							pow(mpion,2)*(1.*pow(mpion,4) + 1.*pow(mrho,4) +
							pow(mpion,2)*(-1.*pow(mrho,2) + 0.5*s) - 0.5*pow(mrho,2)*s) +
							                pow(ma1,2)*(-2.*pow(mpion,4) - 0.5*pow(mrho,2)*s
							+ pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*log(fabs(-1.*pow(ma1,2) + tmin)))/(pow(ma1,2) - 1.*pow(mpion,2)) +
							        (0.25*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-0.5*pow(ma1,6) -
							0.5*pow(mpion,6) + pow(ma1,4)*(0.5*pow(mpion,2) +
							0.5*pow(mrho,2)) +
							                pow(mpion,4)*(0.5*pow(mrho,2) - 1.*s) +
							pow(mpion,2)*(-0.5*pow(mrho,2) + 0.5*s)*s +
							                pow(ma1,2)*(0.5*pow(mpion,4) +
							pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + (-0.5*pow(mrho,2) +
							0.5*s)*s)) +
							             eta1*(1.*pow(mpion,6) + pow(ma1,4)*(1.*pow(mpion,2)
							- 0.5*s) + pow(mpion,2)*(1.5*pow(mrho,2) - 2.*s)*s +
							(-0.5*pow(mrho,2) + 0.5*s)*pow(s,2) +
							                pow(mpion,4)*(-1.*pow(mrho,2) + 1.5*s) +
							pow(ma1,2)*(-2.*pow(mpion,4) + 0.5*pow(mrho,2)*s +
							pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*
							           log(fabs(-1.*pow(ma1,2) + tmin)))/(1.*pow(ma1,2) -
							1.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s) -
							        (0.03125*pow(eta1 - 1.*eta2,2)*(1.*pow(ma1,2) -
							1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s)*
							           (eta1*eta2*(-2.*pow(ma1,8) +
							pow(ma1,6)*(8.*pow(mpion,2) + 4.*pow(mrho,2) - 4.*s) +
							                pow(ma1,4)*(-12.*pow(mpion,4) - 4.*pow(mrho,4) +
							4.*pow(mrho,2)*s + pow(mpion,2)*(-12.*pow(mrho,2) + 8.*s)) +
							                pow(mpion,2)*(-2.*pow(mpion,6) -
							4.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 4.*pow(mrho,4)*s
							- 2.*pow(mrho,2)*pow(s,2) +
							                   pow(mpion,2)*(-8.*pow(mrho,4) +
							8.*pow(mrho,2)*s)) +
							                pow(ma1,2)*(8.*pow(mpion,6) + 2.*pow(mrho,6) +
							pow(mpion,4)*(12.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,4)*s -
							2.*pow(mrho,2)*pow(s,2) + 2.*pow(s,3) +
							                   pow(mpion,2)*(8.*pow(mrho,4) -
							4.*pow(mrho,2)*s - 4.*pow(s,2)))) +
							             pow(eta2,2)*(pow(ma1,8) +
							pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
							                pow(mpion,4)*(pow(mpion,4) +
							2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 1.*pow(mrho,2)*s) +
							                pow(ma1,4)*(6.*pow(mpion,4) - 1.*pow(mrho,4) +
							pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) + pow(mrho,2)*s) +
							                pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mrho,6) -
							5.*pow(mrho,4)*s + 4.*pow(mrho,2)*pow(s,2) - 1.*pow(s,3) +
							pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s) +
							                   pow(mpion,2)*(2.*pow(mrho,4) -
							4.*pow(mrho,2)*s + 2.*pow(s,2)))) +
							             pow(eta1,2)*(pow(ma1,8) + pow(mpion,8) +
							2.*pow(mpion,6)*pow(mrho,2) +
							pow(mpion,2)*pow(mrho,2)*s*(-2.*pow(mrho,2) + 2.*s) +
							                pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) +
							2.*s) + pow(mpion,4)*(3.*pow(mrho,4) - 5.*pow(mrho,2)*s) +
							                pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) +
							pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 3.*pow(mrho,2)*s) +
							                pow(ma1,2)*(-4.*pow(mpion,6) + pow(mrho,4)*s -
							1.*pow(s,3) + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s) +
							                   pow(mpion,2)*(-2.*pow(mrho,4) +
							4.*pow(mrho,2)*s + 2.*pow(s,2)))))*log(fabs(-1.*pow(ma1,2) + tmin)))/
							         (0.25*pow(Gammaa1,2)*pow(ma1,2) + 1.*pow(ma1,4) +
							1.*pow(mpion,4) + 1.*pow(mpion,2)*pow(mrho,2) + 0.25*pow(mrho,4)
							- 1.*pow(mpion,2)*s - 0.5*pow(mrho,2)*s +
							           0.25*pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) -
							1.*pow(mrho,2) + 1.*s)) -
							        (1.*(1.*eta1 - 1.*eta2)*(eta2*(pow(ma1,4)*(0.5*pow(mrho,2)
							- 1.*C4*pow(mrho,4)) +
							                pow(mpion,2)*pow(mrho,2)*(pow(mpion,2)*(0.5 -
							1.*C4*pow(mrho,2)) + (-0.25 + 0.125*delta)*(pow(mrho,2) + s)) +
							                pow(ma1,2)*(-1.*C4*pow(mrho,6) +
							pow(mpion,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4)) +
							0.25*delta*pow(s,2) +
							                   pow(mrho,2)*s*(-0.25 - 0.375*delta - 1.*C4*s) +
							pow(mrho,4)*(0.75 - 0.125*delta + 2.*C4*s))) +
							             eta1*(pow(ma1,4)*(-0.5*pow(mrho,2) +
							1.*C4*pow(mrho,4)) +
							                pow(ma1,2)*(0.5*pow(mrho,4) - 1.*C4*pow(mrho,6) +
							pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) -
							0.25*delta*pow(s,2) +
							                   1.*C4*pow(mrho,2)*pow(s,2)) +
							pow(mrho,2)*(pow(mpion,4)*(-0.5 + 1.*C4*pow(mrho,2)) + s*((0.25 -
							0.125*delta)*pow(mrho,2) + (-0.25 + 0.125*delta)*s) +
							                   pow(mpion,2)*(2.*C4*pow(mrho,4) + (0.5 +
							0.25*delta)*s + pow(mrho,2)*(-1. - 2.*C4*s)))))*log(fabs(-1.*pow(ma1,2) + tmin)))/pow(mrho,2) +
							        0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-1.*pow(mpion,2) + tmin)) +
							        (0.25*(0. + 8.000000000000002*pow(2. -
							1.*delta,2)*pow(mpion,4)*pow(mrho,2) - 5.999999999999999*pow(2. -
							1.*delta,2)*pow(mpion,2)*pow(mrho,2)*s +
							             1.*pow(2. -
							1.*delta,2)*pow(mrho,2)*pow(s,2))*log(fabs(-1.*pow(mpion,2) + tmin)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
							        (0.125*(-2. + delta)*(eta1 - 1.*eta2)*pow(mpion,2)*(0. +
							eta2*pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) +
							             eta1*(2.*pow(mrho,4) - 2.*pow(mrho,2)*s +
							pow(mpion,2)*(-4.*pow(mrho,2) + 4.*s)))*log(fabs(-1.*pow(mpion,2) + tmin)))/(pow(ma1,2) - 1.*pow(mpion,2)) +
							        (2.*(-2. + 1.*delta)*(0. + (-0.25 +
							0.125*delta)*pow(mrho,2)*s + pow(mpion,2)*(-2.*C4*pow(mrho,4) -
							0.5*delta*s + pow(mrho,2)*(1. + 2.*C4*s)))*
							           log(fabs(-1.*pow(mpion,2) + tmin)))/pow(mrho,2) - (0.5*(-2. +
							delta)*(eta1 - 1.*eta2)*pow(mpion,2)*
							           (eta1*(pow(mpion,4)*(-1.*pow(mrho,2) + 1.*s) +
							pow(ma1,2)*(pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) +
							(-0.5*pow(mrho,2) + 0.5*s)*s) +
							                pow(mpion,2)*(-1.*pow(mrho,4) + 2.5*pow(mrho,2)*s
							- 1.5*pow(s,2)) + s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s +
							0.5*pow(s,2))) +
							             eta2*(0.5*pow(mrho,6) + pow(mpion,4)*(1.*pow(mrho,2)
							- 1.*s) - 1.5*pow(mrho,4)*s + 1.5*pow(mrho,2)*pow(s,2) -
							0.5*pow(s,3) +
							                pow(mpion,2)*(1.5*pow(mrho,4) - 3.*pow(mrho,2)*s
							+ 1.5*pow(s,2)) +
							                pow(ma1,2)*(-0.5*pow(mrho,4) + 1.*pow(mrho,2)*s -
							0.5*pow(s,2) + pow(mpion,2)*(-1.*pow(mrho,2) +
							1.*s))))*log(fabs(-1.*pow(mpion,2) + tmin)))/
							         (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4)
							+ 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 2.*pow(mpion,2)*s
							- 2.*pow(mrho,2)*s + pow(s,2) +
							           pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s))
							- 0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-1.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmin)) +
							        (0.5*(-2. + delta)*(eta1 -
							1.*eta2)*(eta2*pow(mpion,6)*(1.*pow(mrho,2) - 1.*s) +
							eta2*pow(ma1,2)*pow(mpion,4)*(-1.*pow(mrho,2) + 1.*s) +
							             eta1*pow(ma1,2)*pow(mpion,2)*(-0.5*pow(mrho,4) +
							pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + 0.5*pow(mrho,2)*s) +
							             eta1*pow(mpion,4)*(0.5*pow(mrho,4) -
							0.5*pow(mrho,2)*s + pow(mpion,2)*(-1.*pow(mrho,2) +
							1.*s)))*log(fabs(-1.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmin)))/
							         (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) -
							2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) +
							        (0.5*(-2. + delta)*(eta1 -
							1.*eta2)*pow(mpion,2)*(eta1*(pow(mpion,2)*(1.*pow(mrho,2) - 1.*s)
							+ (-0.5*pow(mrho,2) + 0.5*s)*s) +
							             eta2*(-0.5*pow(mrho,4) + 1.*pow(mrho,2)*s -
							0.5*pow(s,2) + pow(mpion,2)*(-1.*pow(mrho,2) +
							1.*s)))*log(fabs(-1.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmin)))/
							         (1.*pow(ma1,2) - 1.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s) -
							        (0.25*(0. + 8.000000000000002*pow(2. -
							1.*delta,2)*pow(mpion,4)*pow(mrho,2) + 1.*pow(2. -
							1.*delta,2)*pow(mrho,4)*s +
							             pow(mpion,2)*(C4*(32. - 16.*delta)*pow(mrho,6) +
							delta*(-8. + 4.*delta)*pow(s,2) +
							                pow(mrho,2)*s*(-8. + 24.*delta - 10.*pow(delta,2) +
							32.*C4*s - 16.*C4*delta*s) + pow(mrho,4)*(-16. + 8.*delta - 64.*C4*s +
							32.*C4*delta*s)))*
							           log(fabs(-1.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmin)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
							        0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(4.*pow(ma1,6) +
							pow(Gammaa1,2)*pow(ma1,2)*(-4.*pow(ma1,2) + 4.*pow(mpion,2) -
							2.*s) +
							              pow(ma1,4)*(-12.*pow(mpion,2) + 6.*s) +
							pow(mpion,2)*(-4.*pow(mpion,4) + 4.*pow(mrho,4) +
							2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s) +
							              pow(ma1,2)*(12.*pow(mpion,4) - 2.*pow(mrho,4) -
							8.*pow(mpion,2)*s - 4.*pow(mrho,2)*s + 4.*pow(s,2))) +
							           pow(eta1,2)*(-2.*pow(ma1,6) + 2.*pow(mpion,6) +
							3.*pow(mpion,4)*pow(mrho,2) + pow(ma1,4)*(6.*pow(mpion,2) +
							3.*pow(mrho,2) - 3.*s) - 1.*pow(mpion,4)*s -
							              1.*pow(mpion,2)*pow(mrho,2)*s - 1.*pow(mrho,4)*s +
							pow(mrho,2)*pow(s,2) +
							              pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,2) -
							2.*pow(mpion,2) - 1.*pow(mrho,2) + s) +
							              pow(ma1,2)*(-6.*pow(mpion,4) - 1.*pow(mrho,4) +
							4.*pow(mrho,2)*s - 2.*pow(s,2) + pow(mpion,2)*(-6.*pow(mrho,2) +
							4.*s))) +
							           pow(eta2,2)*(-2.*pow(ma1,6) + 2.*pow(mpion,6) -
							3.*pow(mpion,4)*pow(mrho,2) + pow(ma1,4)*(6.*pow(mpion,2) -
							3.*pow(mrho,2) - 3.*s) - 1.*pow(mpion,4)*s +
							              pow(mpion,2)*pow(mrho,2)*s +
							pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,2) - 2.*pow(mpion,2) +
							pow(mrho,2) + s) +
							              pow(ma1,2)*(-6.*pow(mpion,4) - 1.*pow(mrho,4) +
							2.*pow(mrho,2)*s - 2.*pow(s,2) + pow(mpion,2)*(6.*pow(mrho,2) +
							4.*s))))*
							         log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) -
							4.*pow(ma1,2)*pow(mpion,2) + 4.*pow(mpion,4) +
							2.*pow(ma1,2)*(-1.*pow(mrho,2) + s + tmin) -
							           4.*pow(mpion,2)*(-1.*pow(mrho,2) + s + tmin) +
							pow(-1.*pow(mrho,2) + s + tmin,2))) -
							        (0.5*(1.*eta1 -
							1.*eta2)*(eta2*(pow(Gammaa1,2)*pow(ma1,2)*pow(mrho,2)*(0.5 -
							1.*C4*pow(mrho,2)) + pow(ma1,4)*(-0.5*pow(mrho,2) +
							1.*C4*pow(mrho,4)) +
							                pow(mpion,2)*pow(mrho,2)*(pow(mpion,2)*(-0.5 +
							1.*C4*pow(mrho,2)) + (0.25 - 0.125*delta)*(pow(mrho,2) + s)) +
							                pow(ma1,2)*(1.*C4*pow(mrho,6) +
							pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) -
							0.25*delta*pow(s,2) + pow(mrho,4)*(-0.75 + 0.125*delta - 2.*C4*s) +
							                   pow(mrho,2)*s*(0.25 + 0.375*delta + 1.*C4*s))) + eta1*
							              (pow(Gammaa1,2)*pow(ma1,2)*pow(mrho,2)*(-0.5 +
							1.*C4*pow(mrho,2)) + pow(ma1,4)*(0.5*pow(mrho,2) -
							1.*C4*pow(mrho,4)) +
							                pow(ma1,2)*(-0.5*pow(mrho,4) + 1.*C4*pow(mrho,6)
							+ pow(mpion,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4)) +
							0.25*delta*pow(s,2) -
							                   1.*C4*pow(mrho,2)*pow(s,2)) +
							pow(mrho,2)*(pow(mpion,4)*(0.5 - 1.*C4*pow(mrho,2)) + s*((-0.25 +
							0.125*delta)*pow(mrho,2) + (0.25 - 0.125*delta)*s) +
							                   pow(mpion,2)*(-2.*C4*pow(mrho,4) + (-0.5 -
							0.25*delta)*s + pow(mrho,2)*(1. + 2.*C4*s)))))*
							           log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) -
							2.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmin,2))))/pow(mrho,2) +
							        (0.125*(-2. + delta)*(eta1 -
							1.*eta2)*(eta2*(0.5*pow(Gammaa1,4)*pow(ma1,4) - 0.5*pow(ma1,8) +
							0.5*pow(mpion,8) + 0.5*pow(ma1,4)*pow(mpion,2)*pow(mrho,2) -
							                0.5*pow(mpion,6)*pow(mrho,2) +
							pow(Gammaa1,2)*(pow(ma1,2)*pow(mpion,2)*(1.*pow(mpion,2) +
							1.5*pow(mrho,2) - 2.*s) +
							                   pow(ma1,4)*(-1.*pow(mpion,2) + 0.5*pow(mrho,2)
							- 1.*s)) + pow(ma1,6)*(1.*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) +
							                pow(ma1,2)*pow(mpion,4)*(-1.*pow(mpion,2) -
							0.5*pow(mrho,2) + 1.*s)) +
							             eta1*(pow(ma1,6)*(1.*pow(mpion,2) + 0.5*s) +
							pow(ma1,2)*(3.*pow(mpion,6) + 1.*pow(mpion,2)*pow(mrho,4) -
							0.5*pow(mpion,4)*s) +
							                pow(ma1,4)*(-3.*pow(mpion,4) +
							pow(mpion,2)*(-1.*pow(mrho,2) + 0.5*s) - 0.5*pow(mrho,2)*s) +
							                pow(mpion,4)*(-1.*pow(mpion,4) - 1.*pow(mrho,4) +
							pow(mpion,2)*(1.*pow(mrho,2) - 0.5*s) + 0.5*pow(mrho,2)*s) +
							                pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(mpion,4) +
							pow(ma1,2)*(1.*pow(mpion,2) + 0.5*s) - 0.5*pow(mrho,2)*s +
							pow(mpion,2)*(-1.*pow(mrho,2) + 1.5*s))))*
							           log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) +
							4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
							4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s +
							             pow(s,2) - 4.*pow(mpion,2)*tmin - 2.*pow(mrho,2)*tmin +
							2.*s*tmin + pow(tmin,2) + pow(ma1,2)*(-4.*pow(mpion,2) -
							2.*pow(mrho,2) + 2.*s + 2.*tmin))))/
							         (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) -
							2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) -
							        (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta1*(-1.*pow(mpion,8)
							+ 1.*pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(1.*pow(mpion,2) -
							0.5*s) - 0.5*pow(mpion,6)*s -
							                4.*pow(mpion,4)*pow(mrho,2)*s -
							1.5*pow(mpion,2)*pow(mrho,4)*s + 3.5*pow(mpion,4)*pow(s,2) +
							4.*pow(mpion,2)*pow(mrho,2)*pow(s,2) +
							                0.5*pow(mrho,4)*pow(s,2) -
							2.5*pow(mpion,2)*pow(s,3) - 1.*pow(mrho,2)*pow(s,3) +
							0.5*pow(s,4) +
							                pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(mpion,4) +
							pow(ma1,2)*(1.*pow(mpion,2) - 0.5*s) - 0.5*pow(mpion,2)*s +
							0.5*pow(s,2)) +
							                pow(ma1,4)*(-3.*pow(mpion,4) + (1.*pow(mrho,2) -
							0.5*s)*s + pow(mpion,2)*(-2.*pow(mrho,2) + 2.5*s)) +
							                pow(ma1,2)*(3.*pow(mpion,6) +
							pow(mpion,4)*(2.*pow(mrho,2) - 1.5*s) - 0.5*pow(mrho,4)*s +
							0.5*pow(s,3) +
							                   pow(mpion,2)*(1.*pow(mrho,4) -
							1.*pow(mrho,2)*s - 1.*pow(s,2)))) +
							             eta2*(0.5*pow(Gammaa1,4)*pow(ma1,4) - 0.5*pow(ma1,8)
							+ 0.5*pow(mpion,8) + pow(ma1,6)*(1.*pow(mpion,2) +
							1.*pow(mrho,2) - 0.5*s) + 0.5*pow(mpion,6)*s +
							                pow(ma1,4)*(-0.5*pow(mrho,4) + (-0.5*pow(mpion,2)
							+ 0.5*s)*s) + pow(mpion,4)*(-0.5*pow(mrho,4) + 2.*pow(mrho,2)*s -
							1.5*pow(s,2)) +
							                pow(mpion,2)*s*(0.5*pow(mrho,4) -
							1.*pow(mrho,2)*s + 0.5*pow(s,2)) +
							                pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,4) +
							0.5*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 1.5*s) -
							1.*pow(mrho,2)*s + 0.5*pow(s,2) +
							                   pow(ma1,2)*(-1.*pow(mpion,2) - 1.*pow(mrho,2)
							+ 1.5*s)) +
							                pow(ma1,2)*(-1.*pow(mpion,6) +
							pow(mpion,4)*(-1.*pow(mrho,2) + 0.5*s) +
							pow(mpion,2)*(-1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.*pow(s,2)) +
							                   s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s +
							0.5*pow(s,2)))))*
							           log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) +
							4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
							4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s +
							             pow(s,2) - 4.*pow(mpion,2)*tmin - 2.*pow(mrho,2)*tmin +
							2.*s*tmin + pow(tmin,2) + pow(ma1,2)*(-4.*pow(mpion,2) -
							2.*pow(mrho,2) + 2.*s + 2.*tmin))))/
							         (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4)
							+ 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 2.*pow(mpion,2)*s
							- 2.*pow(mrho,2)*s + pow(s,2) +
							           pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s)) -
							        (0.0625*pow(eta1 -
							1.*eta2,2)*(pow(eta2,2)*(-1.*pow(ma1,10) +
							pow(ma1,8)*(5.*pow(mpion,2) + 2.5*pow(mrho,2) - 2.5*s) +
							                pow(Gammaa1,4)*pow(ma1,4)*(1.*pow(ma1,2) -
							1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s) +
							                pow(ma1,4)*(10.*pow(mpion,6) - 2.5*pow(mrho,6) +
							pow(mpion,4)*(15.*pow(mrho,2) - 9.*s) + 6.*pow(mrho,4)*s -
							4.5*pow(mrho,2)*pow(s,2) + 1.*pow(s,3)) +
							                pow(ma1,6)*(-10.*pow(mpion,4) + (1.*pow(mrho,2) -
							1.*s)*s + pow(mpion,2)*(-10.*pow(mrho,2) + 8.*s)) +
							                pow(mpion,4)*(1.*pow(mpion,6) + 0.5*pow(mrho,6) +
							pow(mpion,4)*(2.5*pow(mrho,2) - 0.5*s) - 1.*pow(mrho,4)*s +
							0.5*pow(mrho,2)*pow(s,2) +
							                   pow(mpion,2)*(2.*pow(mrho,4) -
							2.*pow(mrho,2)*s)) +
							                pow(Gammaa1,2)*pow(ma1,2)*(4.*pow(ma1,6) -
							4.*pow(mpion,6) - 0.5*pow(mrho,6) + 1.5*pow(mrho,4)*s -
							1.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) +
							                   pow(mpion,4)*(-6.*pow(mrho,2) + 6.*s) +
							pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
							                   pow(mpion,2)*(-3.*pow(mrho,4) +
							6.*pow(mrho,2)*s - 3.*pow(s,2)) +
							                   pow(ma1,2)*(12.*pow(mpion,4) + 3.*pow(mrho,4)
							+ pow(mpion,2)*(12.*pow(mrho,2) - 12.*s) - 6.*pow(mrho,2)*s +
							3.*pow(s,2))) +
							                pow(ma1,2)*(-5.*pow(mpion,8) + 1.*pow(mrho,8) -
							3.5*pow(mrho,6)*s + 4.5*pow(mrho,4)*pow(s,2) -
							2.5*pow(mrho,2)*pow(s,3) + 0.5*pow(s,4) +
							                   pow(mpion,6)*(-10.*pow(mrho,2) + 4.*s) +
							pow(mpion,4)*(-2.*pow(mrho,4) + 1.*pow(mrho,2)*s + 1.*pow(s,2)) +
							                   pow(mpion,2)*(3.*pow(mrho,6) -
							8.*pow(mrho,4)*s + 7.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3)))) +
							             pow(eta1,2)*(-1.*pow(ma1,10) +
							pow(ma1,8)*(5.*pow(mpion,2) + 2.5*pow(mrho,2) - 2.5*s) +
							                pow(Gammaa1,4)*pow(ma1,4)*(1.*pow(ma1,2) -
							1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s) +
							                pow(ma1,6)*(-10.*pow(mpion,4) - 2.*pow(mrho,4) +
							5.*pow(mrho,2)*s - 1.*pow(s,2) + pow(mpion,2)*(-10.*pow(mrho,2)
							+ 8.*s)) +
							                pow(ma1,4)*(10.*pow(mpion,6) + 0.5*pow(mrho,6) +
							pow(mpion,4)*(15.*pow(mrho,2) - 9.*s) - 3.*pow(mrho,4)*s +
							1.5*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) +
							                   pow(mpion,2)*(6.*pow(mrho,4) -
							12.*pow(mrho,2)*s)) +
							                pow(Gammaa1,2)*pow(ma1,2)*(4.*pow(ma1,6) -
							4.*pow(mpion,6) - 0.5*pow(mrho,6) + 1.5*pow(mrho,4)*s -
							1.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) +
							                   pow(mpion,4)*(-6.*pow(mrho,2) + 6.*s) +
							pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
							                   pow(mpion,2)*(-3.*pow(mrho,4) +
							6.*pow(mrho,2)*s - 3.*pow(s,2)) +
							                   pow(ma1,2)*(12.*pow(mpion,4) + 3.*pow(mrho,4)
							+ pow(mpion,2)*(12.*pow(mrho,2) - 12.*s) - 6.*pow(mrho,2)*s +
							3.*pow(s,2))) +
							                pow(mpion,2)*(1.*pow(mpion,8) +
							pow(mpion,6)*(2.5*pow(mrho,2) - 0.5*s) +
							pow(mpion,4)*(4.*pow(mrho,4) - 6.*pow(mrho,2)*s) +
							                   pow(mrho,2)*s*(-1.*pow(mrho,4) +
							2.*pow(mrho,2)*s - 1.*pow(s,2)) + pow(mpion,2)*(1.5*pow(mrho,6)
							- 6.*pow(mrho,4)*s + 4.5*pow(mrho,2)*pow(s,2))) +
							                pow(ma1,2)*(-5.*pow(mpion,8) +
							pow(mpion,6)*(-10.*pow(mrho,2) + 4.*s) +
							pow(mpion,4)*(-8.*pow(mrho,4) + 13.*pow(mrho,2)*s + 1.*pow(s,2)) +
							                   pow(mpion,2)*(-1.*pow(mrho,6) +
							6.*pow(mrho,4)*s - 3.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3)) +
							                   s*(0.5*pow(mrho,6) - 0.5*pow(mrho,4)*s -
							0.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3)))) +
							             eta1*eta2*(2.*pow(ma1,10) +
							pow(Gammaa1,4)*pow(ma1,4)*(-2.*pow(ma1,2) + 2.*pow(mpion,2) +
							1.*pow(mrho,2) - 1.*s) +
							                pow(ma1,8)*(-10.*pow(mpion,2) - 5.*pow(mrho,2) +
							5.*s) +
							                pow(ma1,6)*(20.*pow(mpion,4) + 6.*pow(mrho,4) +
							pow(mpion,2)*(20.*pow(mrho,2) - 16.*s) - 8.*pow(mrho,2)*s +
							2.*pow(s,2)) +
							                pow(ma1,4)*(-20.*pow(mpion,6) - 4.*pow(mrho,6) +
							6.*pow(mrho,4)*s - 2.*pow(s,3) + pow(mpion,4)*(-30.*pow(mrho,2)
							+ 18.*s) +
							                   pow(mpion,2)*(-18.*pow(mrho,4) +
							18.*pow(mrho,2)*s)) +
							                pow(mpion,2)*(-2.*pow(mpion,8) - 1.*pow(mrho,8) +
							3.*pow(mrho,6)*s - 3.*pow(mrho,4)*pow(s,2) +
							1.*pow(mrho,2)*pow(s,3) +
							                   pow(mpion,6)*(-5.*pow(mrho,2) + 1.*s) +
							pow(mpion,4)*(-10.*pow(mrho,4) + 10.*pow(mrho,2)*s) +
							                   pow(mpion,2)*(-6.*pow(mrho,6) +
							12.*pow(mrho,4)*s - 6.*pow(mrho,2)*pow(s,2))) +
							                pow(ma1,2)*(10.*pow(mpion,8) + 1.*pow(mrho,8) +
							pow(mpion,6)*(20.*pow(mrho,2) - 8.*s) - 2.*pow(mrho,6)*s +
							2.*pow(mrho,2)*pow(s,3) - 1.*pow(s,4) +
							                   pow(mpion,4)*(22.*pow(mrho,4) -
							20.*pow(mrho,2)*s - 2.*pow(s,2)) + pow(mpion,2)*(8.*pow(mrho,6)
							- 12.*pow(mrho,4)*s + 4.*pow(s,3))) +
							                pow(Gammaa1,2)*pow(ma1,2)*(-8.*pow(ma1,6) +
							8.*pow(mpion,6) + 1.*pow(mrho,6) + pow(mpion,4)*(12.*pow(mrho,2)
							- 12.*s) +
							                   pow(ma1,4)*(24.*pow(mpion,2) + 12.*pow(mrho,2)
							- 12.*s) - 3.*pow(mrho,4)*s + 3.*pow(mrho,2)*pow(s,2) -
							1.*pow(s,3) +
							                   pow(mpion,2)*(6.*pow(mrho,4) -
							12.*pow(mrho,2)*s + 6.*pow(s,2)) +
							                   pow(ma1,2)*(-24.*pow(mpion,4) - 6.*pow(mrho,4)
							+ 12.*pow(mrho,2)*s - 6.*pow(s,2) +
							pow(mpion,2)*(-24.*pow(mrho,2) + 24.*s)))))*
							           log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) +
							4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
							4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s +
							             pow(s,2) - 4.*pow(mpion,2)*tmin - 2.*pow(mrho,2)*tmin +
							2.*s*tmin + pow(tmin,2) + pow(ma1,2)*(-4.*pow(mpion,2) -
							2.*pow(mrho,2) + 2.*s + 2.*tmin))))/
							         (pow(Gammaa1,2)*pow(ma1,2) + 4.*pow(ma1,4) +
							4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
							4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s +
							           pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) -
							4.*pow(mrho,2) + 4.*s))))/(16.*Pi*s*(-4*pow(mpion,2) + s)) +
							  (pow(Const,2)*pow(ghat,4)*(0. + (0.03125*pow(eta1 -
							1.*eta2,2)*(eta1*eta2*
							             (-2.*pow(ma1,8) - 2.*pow(mpion,8) +
							2.*pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) +
							               pow(ma1,2)*pow(mpion,2)*(8.*pow(mpion,4) -
							8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
							               pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) +
							8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
							            pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) -
							2.*pow(mpion,6)*pow(mrho,2) + 1.*pow(mpion,4)*pow(mrho,4) +
							               pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) +
							2.*s) +
							               pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) +
							pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
							2.*pow(s,2)) +
							               pow(ma1,2)*(-4.*pow(mpion,6) -
							2.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(6.*pow(mrho,2) +
							2.*s))) +
							            pow(eta1,2)*(1.*pow(ma1,8) +
							pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
							               pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) +
							pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
							2.*pow(s,2)) +
							               pow(ma1,2)*(-4.*pow(mpion,6) +
							2.*pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(2.*pow(mrho,2) -
							2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
							               pow(mpion,2)*(1.*pow(mpion,6) +
							2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s
							+ pow(mpion,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s)))))/
							        (1.*pow(ma1,2) - 1.*tmax) + (1.*pow(-2. +
							delta,2)*pow(mpion,2)*(1.*pow(mpion,2) -
							0.25*pow(mrho,2)))/(1.*pow(mpion,2) - 1.*tmax) +
							       (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) -
							0.25*pow(mrho,2)))/(1.*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s -
							1.*tmax) -
							       (0.5*pow(-2. + delta,2)*pow(mpion,2)*tmax)/pow(mrho,2) -
							0.25*(-2. + delta)*(eta1 - 1.*eta2)*
							        (-0.5*eta2*pow(ma1,2) + 1.*eta1*pow(mpion,2) +
							0.5*eta2*pow(mrho,2) + 0.5*eta1*s - 1.*eta2*s)*tmax +
							       (0.25*(pow(mpion,2)*(12. + 1.*pow(delta,2) -
							16.*C4*pow(mrho,2) + delta*(-8. + 8.*C4*pow(mrho,2))) +
							            (-4. - 3.*pow(delta,2) - 16.*C4*pow(mrho,2) + delta*(8.
							+ 8.*C4*pow(mrho,2)))*s)*tmax)/pow(mrho,2) -
							       0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(pow(ma1,2) - 1.*s)
							+ eta1*(-2.*pow(mpion,2) + s))*tmax -
							       0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-1.*pow(ma1,2) -
							2.*pow(mpion,2) - 1.*s) + eta1*(2.*pow(mpion,2) + s))*tmax +
							       0.25*(-2. + delta)*(eta1 - 1.*eta2)*(eta1*(1.*pow(mpion,2) -
							0.5*s) + eta2*(-0.5*pow(ma1,2) - 1.*pow(mpion,2) - 0.5*pow(mrho,2)
							+ 1.*s))*tmax +
							       (0.25*(-2. + 1.*delta)*(-8.*C4*pow(mrho,4) +
							pow(mpion,2)*(2. + 1.*delta - 8.*C4*pow(mrho,2)) + (-2. -
							3.*delta)*s + pow(mrho,2)*(2. + 1.*delta + 16.*C4*s))*tmax)/
							        pow(mrho,2) + (0.25*(32*pow(C4,2)*pow(mrho,8) +
							2*pow(delta,2)*pow(s,2) + 8*C4*pow(mrho,6)*(-6 + delta - 8*C4*s) +
							2*delta*pow(mrho,2)*s*(-6 + delta - 8*C4*s) +
							            pow(mrho,4)*(12 - pow(delta,2) + 8*C4*(6 + delta)*s +
							32*pow(C4,2)*pow(s,2)))*tmax)/pow(mrho,4) -
							       (1.*(1.*eta1 - 1.*eta2)*(eta2*(0.75*pow(mrho,4) -
							0.125*delta*pow(mrho,4) - 1.*C4*pow(mrho,6) +
							pow(ma1,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) +
							               pow(mpion,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4))
							- 0.25*pow(mrho,2)*s - 0.375*delta*pow(mrho,2)*s +
							2.*C4*pow(mrho,4)*s + 0.25*delta*pow(s,2) -
							               1.*C4*pow(mrho,2)*pow(s,2)) +
							eta1*(0.5*pow(mrho,4) - 1.*C4*pow(mrho,6) +
							pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) +
							               pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4))
							- 0.25*delta*pow(s,2) +
							1.*C4*pow(mrho,2)*pow(s,2)))*tmax)/pow(mrho,2) +
							       0.0625*pow(eta1 -
							1.*eta2,2)*(pow(eta2,2)*(pow(Gammaa1,2)*pow(ma1,2) -
							1.*pow(ma1,4) - 2.*pow(mpion,4) - 2.*pow(mpion,2)*pow(mrho,2) +
							2.*pow(mrho,4) +
							             pow(ma1,2)*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s) -
							3.*pow(mrho,2)*s + pow(s,2)) +
							          pow(eta1,2)*(pow(Gammaa1,2)*pow(ma1,2) -
							1.*pow(ma1,4) - 2.*pow(mpion,4) - 2.*pow(mpion,2)*pow(mrho,2) +
							             pow(ma1,2)*(2.*pow(mpion,2) + pow(mrho,2) - 1.*s) +
							pow(mrho,2)*s + pow(s,2)) +
							          eta1*eta2*(-2.*pow(Gammaa1,2)*pow(ma1,2) +
							2.*pow(ma1,4) + 4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) +
							2.*pow(mrho,4) - 2.*pow(s,2) +
							             pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) +
							2.*s)))*tmax +
							       0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-6.*pow(ma1,4) -
							12.*pow(mpion,4) + 2.*pow(mrho,4) + pow(ma1,2)*(16.*pow(mpion,2)
							- 8.*s) + 8.*pow(mpion,2)*s +
							             4.*pow(mrho,2)*s - 4.*pow(s,2)) +
							pow(eta1,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) +
							pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
							             2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) -
							4.*pow(mrho,2) + 4.*s)) +
							          pow(eta2,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) +
							pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) -
							2.*pow(mrho,2)*s + 2.*pow(s,2) +
							             pow(ma1,2)*(-8.*pow(mpion,2) + 4.*pow(mrho,2) +
							4.*s)))*tmax - (0.125*(-2. + 1.*delta)*(2. + 1.*delta -
							8.*C4*pow(mrho,2))*pow(tmax,2))/pow(mrho,2) -
							       0.5*pow(1.*eta1 - 1.*eta2,2)*(-0.5 +
							1.*C4*pow(mrho,2))*pow(tmax,2) -
							       (1.*(0.5 - 0.125*pow(delta,2) - 2.*C4*pow(mrho,2) +
							1.*C4*delta*pow(mrho,2))*pow(tmax,2))/pow(mrho,2) +
							       0.0625*pow(1.*eta1 - 1.*eta2,4)*(1.*pow(mpion,2) +
							0.5*pow(mrho,2) - 0.5*s)*pow(tmax,2) +
							       0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) +
							2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) + eta1*(pow(ma1,2) -
							2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*
							        pow(tmax,2) + 0.010416666666666666*pow(eta1 -
							1.*eta2,4)*pow(tmax,3) - 0.020833333333333332*pow(1.*eta1 -
							1.*eta2,4)*pow(tmax,3) +
							       0.03125*pow(eta1 -
							1.*eta2,2)*(eta1*eta2*(2.*pow(Gammaa1,2)*pow(ma1,2) -
							6.*pow(ma1,4) - 4.*pow(mpion,4) + 2.*pow(mrho,4) +
							pow(ma1,2)*(8.*pow(mpion,2) - 8.*s) +
							             4.*pow(mrho,2)*s - 4.*pow(s,2)) +
							pow(eta1,2)*(-1.*pow(Gammaa1,2)*pow(ma1,2) + 3.*pow(ma1,4) +
							2.*pow(mpion,4) + 2.*pow(mpion,2)*pow(mrho,2) +
							             pow(mrho,4) - 4.*pow(mrho,2)*s + 2.*pow(s,2) +
							pow(ma1,2)*(-4.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
							          pow(eta2,2)*(-1.*pow(Gammaa1,2)*pow(ma1,2) +
							3.*pow(ma1,4) + 2.*pow(mpion,4) - 2.*pow(mpion,2)*pow(mrho,2) +
							pow(mrho,4) - 2.*pow(mrho,2)*s +
							             2.*pow(s,2) + pow(ma1,2)*(-4.*pow(mpion,2) +
							4.*pow(mrho,2) + 4.*s)))*(-1.*pow(mrho,2) + s + tmax) -
							       0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) -
							1.*pow(mrho,2) - 1.*s) + eta1*(pow(ma1,2) - 1.*pow(mrho,2) +
							s))*pow(-1.*pow(mrho,2) + s + tmax,2) +
							       0.010416666666666666*pow(eta1 -
							1.*eta2,4)*pow(-1.*pow(mrho,2) + s + tmax,3) +
							       0.25*(eta1 - 1.*eta2)*(1.*eta1 - 1.*eta2)*(-1. +
							2.*C4*pow(mrho,2))*pow(pow(ma1,2) - 2.*pow(mpion,2) -
							1.*pow(mrho,2) + s + tmax,2) -
							       (2.*(1.*eta1 - 1.*eta2)*(eta2*(0.375*pow(mrho,4) -
							0.0625*delta*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
							pow(ma1,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) +
							               pow(mpion,2)*(-0.5*pow(mrho,2) +
							1.*C4*pow(mrho,4)) - 0.125*pow(mrho,2)*s -
							0.1875*delta*pow(mrho,2)*s + 1.*C4*pow(mrho,4)*s +
							0.125*delta*pow(s,2) -
							               0.5*C4*pow(mrho,2)*pow(s,2)) +
							eta1*(0.25*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
							pow(mpion,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) +
							               pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4))
							- 0.125*delta*pow(s,2) + 0.5*C4*pow(mrho,2)*pow(s,2)))*
							          (1.*pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) +
							1.*s + 1.*tmax))/pow(mrho,2) +
							       (2.*(1.*eta1 - 1.*eta2)*Gammaa1*ma1*(eta2*(0.375*pow(mrho,4) -
							0.0625*delta*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
							pow(ma1,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) +
							               pow(mpion,2)*(-0.5*pow(mrho,2) +
							1.*C4*pow(mrho,4)) - 0.125*pow(mrho,2)*s -
							0.1875*delta*pow(mrho,2)*s + 1.*C4*pow(mrho,4)*s +
							0.125*delta*pow(s,2) -
							               0.5*C4*pow(mrho,2)*pow(s,2)) +
							eta1*(0.25*pow(mrho,4) - 0.5*C4*pow(mrho,6) +
							pow(mpion,2)*(0.5*pow(mrho,2) - 1.*C4*pow(mrho,4)) +
							               pow(ma1,2)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4))
							- 0.125*delta*pow(s,2) + 0.5*C4*pow(mrho,2)*pow(s,2)))*
							          atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) +
							s + tmax)/(Gammaa1*ma1)))/pow(mrho,2) +
							       (0.25*(-2. + delta)*(eta1 -
							1.*eta2)*Gammaa1*ma1*(eta2*(-1.*pow(ma1,6) +
							pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(ma1,2) + 0.5*pow(mrho,2) -
							1.*s) +
							               pow(ma1,4)*(2.*pow(mpion,2) + 0.5*pow(mrho,2) -
							1.*s) + pow(mpion,4)*(-1.5*pow(mrho,2) + 1.*s) +
							               pow(ma1,2)*pow(mpion,2)*(-1.*pow(mpion,2) -
							1.*pow(mrho,2) + 2.*s)) +
							            eta1*(pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,2) +
							0.5*s) + pow(ma1,4)*(1.*pow(mpion,2) + 0.5*s) +
							pow(ma1,2)*(-2.*pow(mpion,4) - 1.*pow(mpion,2)*s) +
							               pow(mpion,2)*(1.*pow(mpion,4) - 1.*pow(mrho,4) +
							pow(mpion,2)*(2.*pow(mrho,2) - 1.5*s) + 1.*pow(mrho,2)*s)))*
							          atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) +
							s + tmax)/(Gammaa1*ma1)))/
							        (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) -
							2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) -
							       (0.25*(-2. + delta)*(eta1 -
							1.*eta2)*Gammaa1*ma1*(eta2*(-1.*pow(ma1,6) -
							2.*pow(mpion,4)*pow(mrho,2) - 1.*pow(mpion,2)*pow(mrho,4) +
							               pow(ma1,4)*(2.*pow(mpion,2) + 2.*pow(mrho,2) -
							1.5*s) + 2.5*pow(mpion,4)*s + 3.*pow(mpion,2)*pow(mrho,2)*s +
							0.5*pow(mrho,4)*s -
							               2.*pow(mpion,2)*pow(s,2) -
							1.*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) +
							pow(Gammaa1,2)*(-1.*pow(ma1,4) + 0.5*pow(ma1,2)*s) +
							               pow(ma1,2)*(-1.*pow(mpion,4) - 1.*pow(mrho,4) +
							1.*pow(mrho,2)*s + pow(mpion,2)*(-2.*pow(mrho,2) + 1.*s))) +
							            eta1*(1.*pow(mpion,6) + 4.*pow(mpion,4)*pow(mrho,2) +
							1.*pow(mpion,2)*pow(mrho,4) +
							pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,2) - 0.5*s) +
							               pow(ma1,4)*(1.*pow(mpion,2) - 0.5*s) -
							4.5*pow(mpion,4)*s - 4.*pow(mpion,2)*pow(mrho,2)*s -
							0.5*pow(mrho,4)*s + 3.*pow(mpion,2)*pow(s,2) +
							               1.*pow(mrho,2)*pow(s,2) - 0.5*pow(s,3) +
							pow(ma1,2)*(-2.*pow(mpion,4) + (1.*pow(mrho,2) - 1.*s)*s +
							pow(mpion,2)*(-2.*pow(mrho,2) + 3.*s))))*
							          atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) +
							s + tmax)/(Gammaa1*ma1)))/
							        (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4)
							+ 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 2.*pow(mpion,2)*s
							- 2.*pow(mrho,2)*s + pow(s,2) +
							          pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s))
							+ (0.03125*pow(eta1 - 1.*eta2,2)*
							          (pow(eta2,2)*(pow(Gammaa1,4)*pow(ma1,4) + pow(ma1,8)
							+ pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) +
							pow(mpion,4)*pow(mrho,4) +
							               pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) +
							2.*s) +
							               pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) +
							pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
							2.*pow(s,2)) +
							               pow(ma1,2)*(-4.*pow(mpion,6) -
							2.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(6.*pow(mrho,2) +
							2.*s)) +
							               pow(Gammaa1,2)*pow(ma1,2)*(-6.*pow(ma1,4) -
							6.*pow(mpion,4) - 1.*pow(mrho,4) + pow(ma1,2)*(12.*pow(mpion,2)
							- 6.*pow(mrho,2) - 6.*s) +
							                  2.*pow(mrho,2)*s - 2.*pow(s,2) +
							pow(mpion,2)*(6.*pow(mrho,2) + 4.*s))) +
							            eta1*eta2*(-2.*pow(Gammaa1,4)*pow(ma1,4) -
							2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) +
							pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) +
							               pow(ma1,2)*pow(mpion,2)*(8.*pow(mpion,4) -
							8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
							               pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) +
							8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2)) +
							               pow(Gammaa1,2)*pow(ma1,2)*(12.*pow(ma1,4) +
							12.*pow(mpion,4) - 2.*pow(mrho,4) - 8.*pow(mpion,2)*s -
							4.*pow(mrho,2)*s + 4.*pow(s,2) +
							                  pow(ma1,2)*(-24.*pow(mpion,2) + 12.*s))) +
							pow(eta1,2)*
							             (pow(Gammaa1,4)*pow(ma1,4) + pow(ma1,8) +
							pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
							               pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) +
							pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
							2.*pow(s,2)) +
							               pow(ma1,2)*(-4.*pow(mpion,6) +
							2.*pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(2.*pow(mrho,2) -
							2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
							               pow(Gammaa1,2)*pow(ma1,2)*(-6.*pow(ma1,4) -
							6.*pow(mpion,4) - 1.*pow(mrho,4) + pow(ma1,2)*(12.*pow(mpion,2)
							+ 6.*pow(mrho,2) - 6.*s) +
							                  4.*pow(mrho,2)*s - 2.*pow(s,2) +
							pow(mpion,2)*(-6.*pow(mrho,2) + 4.*s)) +
							               pow(mpion,2)*(pow(mpion,6) +
							2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s
							+ pow(mpion,2)*(pow(mrho,4) - 2.*pow(mrho,2)*s))))*
							          atan((pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) +
							s + tmax)/(Gammaa1*ma1)))/(Gammaa1*ma1) -
							       (0.0625*pow(eta1 -
							1.*eta2,2)*Gammaa1*ma1*(eta1*eta2*(-2.*pow(Gammaa1,4)*pow(ma1,4) +
							14.*pow(ma1,8) + 14.*pow(mpion,8) + 28.*pow(mpion,6)*pow(mrho,2) +
							               20.*pow(mpion,4)*pow(mrho,4) +
							10.*pow(mpion,2)*pow(mrho,6) + 2.*pow(mrho,8) -
							16.*pow(mpion,6)*s - 16.*pow(mpion,4)*pow(mrho,2)*s -
							               12.*pow(mpion,2)*pow(mrho,4)*s - 4.*pow(mrho,6)*s
							- 4.*pow(mpion,4)*pow(s,2) -
							6.*pow(mpion,2)*pow(mrho,2)*pow(s,2) + 8.*pow(mpion,2)*pow(s,3) +
							               4.*pow(mrho,2)*pow(s,3) - 2.*pow(s,4) +
							pow(ma1,6)*(-56.*pow(mpion,2) - 28.*pow(mrho,2) + 28.*s) +
							               pow(ma1,4)*(84.*pow(mpion,4) + 24.*pow(mrho,4) +
							pow(mpion,2)*(84.*pow(mrho,2) - 72.*s) - 36.*pow(mrho,2)*s +
							12.*pow(s,2)) +
							               pow(Gammaa1,2)*pow(ma1,2)*(-4.*pow(ma1,4) -
							4.*pow(mpion,4) + pow(ma1,2)*(8.*pow(mpion,2) + 4.*pow(mrho,2) -
							4.*s) + (4.*pow(mrho,2) - 4.*s)*s +
							                  pow(mpion,2)*(-4.*pow(mrho,2) + 8.*s)) +
							pow(ma1,2)*
							                (-56.*pow(mpion,6) - 10.*pow(mrho,6) +
							18.*pow(mrho,4)*s - 6.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3) +
							pow(mpion,4)*(-84.*pow(mrho,2) + 60.*s) +
							                  pow(mpion,2)*(-48.*pow(mrho,4) +
							60.*pow(mrho,2)*s - 12.*pow(s,2)))) +
							            pow(eta1,2)*(1.*pow(Gammaa1,4)*pow(ma1,4) -
							7.*pow(ma1,8) - 7.*pow(mpion,8) - 14.*pow(mpion,6)*pow(mrho,2) -
							7.*pow(mpion,4)*pow(mrho,4) -
							               2.*pow(mpion,2)*pow(mrho,6) +
							pow(ma1,6)*(28.*pow(mpion,2) + 14.*pow(mrho,2) - 14.*s) +
							8.*pow(mpion,6)*s + 11.*pow(mpion,4)*pow(mrho,2)*s +
							               6.*pow(mpion,2)*pow(mrho,4)*s + 1.*pow(mrho,6)*s +
							2.*pow(mpion,4)*pow(s,2) - 1.*pow(mrho,4)*pow(s,2) -
							4.*pow(mpion,2)*pow(s,3) -
							               1.*pow(mrho,2)*pow(s,3) + 1.*pow(s,4) +
							pow(Gammaa1,2)*pow(ma1,2)*
							                (2.*pow(ma1,4) + 2.*pow(mpion,4) + 1.*pow(mrho,4)
							+ pow(mpion,2)*(2.*pow(mrho,2) - 4.*s) - 1.*pow(mrho,2)*s +
							2.*pow(s,2) +
							                  pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) +
							2.*s)) +
							               pow(ma1,4)*(-42.*pow(mpion,4) - 9.*pow(mrho,4) +
							21.*pow(mrho,2)*s - 6.*pow(s,2) + pow(mpion,2)*(-42.*pow(mrho,2)
							+ 36.*s)) +
							               pow(ma1,2)*(28.*pow(mpion,6) + 2.*pow(mrho,6) +
							pow(mpion,4)*(42.*pow(mrho,2) - 30.*s) - 9.*pow(mrho,4)*s +
							6.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) +
							                  pow(mpion,2)*(18.*pow(mrho,4) -
							36.*pow(mrho,2)*s + 6.*pow(s,2)))) +
							            pow(eta2,2)*(1.*pow(Gammaa1,4)*pow(ma1,4) -
							7.*pow(ma1,8) - 7.*pow(mpion,8) - 14.*pow(mpion,6)*pow(mrho,2) -
							1.*pow(mpion,4)*pow(mrho,4) +
							               6.*pow(mpion,2)*pow(mrho,6) + 2.*pow(mrho,8) +
							pow(ma1,6)*(28.*pow(mpion,2) + 14.*pow(mrho,2) - 14.*s) +
							8.*pow(mpion,6)*s -
							               1.*pow(mpion,4)*pow(mrho,2)*s -
							16.*pow(mpion,2)*pow(mrho,4)*s - 7.*pow(mrho,6)*s +
							2.*pow(mpion,4)*pow(s,2) +
							               14.*pow(mpion,2)*pow(mrho,2)*pow(s,2) +
							9.*pow(mrho,4)*pow(s,2) - 4.*pow(mpion,2)*pow(s,3) -
							5.*pow(mrho,2)*pow(s,3) + 1.*pow(s,4) +
							               pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,4) +
							2.*pow(mpion,4) + 3.*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2)
							- 4.*s) - 5.*pow(mrho,2)*s + 2.*pow(s,2) +
							                  pow(ma1,2)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) +
							2.*s)) +
							               pow(ma1,4)*(-42.*pow(mpion,4) - 3.*pow(mrho,4) +
							9.*pow(mrho,2)*s - 6.*pow(s,2) + pow(mpion,2)*(-42.*pow(mrho,2)
							+ 36.*s)) +
							               pow(ma1,2)*(28.*pow(mpion,6) - 4.*pow(mrho,6) +
							pow(mpion,4)*(42.*pow(mrho,2) - 30.*s) + 9.*pow(mrho,4)*s -
							6.*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) +
							                  pow(mpion,2)*(6.*pow(mrho,4) -
							12.*pow(mrho,2)*s + 6.*pow(s,2)))))*atan((pow(ma1,2) -
							2.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmax)/(Gammaa1*ma1)))/
							        (pow(Gammaa1,2)*pow(ma1,2) + 4.*pow(ma1,4) +
							4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
							4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) +
							          pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s))
							+ 0.0625*pow(eta1 - 1.*eta2,2)*
							        (eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2)
							- 6.*s) + pow(mpion,2)*(4.*pow(mpion,4) - 4.*pow(mrho,4) -
							2.*pow(mpion,2)*s + 2.*pow(mrho,2)*s) +
							             pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) +
							8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
							          pow(eta1,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) +
							pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(pow(mrho,2) - 1.*s)*s
							+ pow(mpion,4)*(-3.*pow(mrho,2) + s) +
							             pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
							             pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) +
							pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
							2.*pow(s,2))) +
							          pow(eta2,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) -
							1.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(3.*pow(mrho,2) + s) +
							             pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
							             pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) +
							pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s +
							2.*pow(s,2))))*log(fabs(-1.*pow(ma1,2) + tmax)) -
							       (0.25*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-0.5*pow(ma1,6) -
							0.5*pow(mpion,6) + 0.5*pow(mpion,4)*pow(mrho,2) +
							               pow(ma1,4)*(0.5*pow(mpion,2) + 0.5*pow(mrho,2) -
							1.*s) + pow(ma1,2)*pow(mpion,2)*(0.5*pow(mpion,2) +
							1.*pow(mrho,2) - 1.*s)) +
							            eta1*(pow(ma1,4)*(1.*pow(mpion,2) + 0.5*s) +
							pow(mpion,2)*(1.*pow(mpion,4) + 1.*pow(mrho,4) +
							pow(mpion,2)*(-1.*pow(mrho,2) + 0.5*s) - 0.5*pow(mrho,2)*s) +
							               pow(ma1,2)*(-2.*pow(mpion,4) - 0.5*pow(mrho,2)*s +
							pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*log(fabs(-1.*pow(ma1,2) + tmax)))/(pow(ma1,2) - 1.*pow(mpion,2)) +
							       (0.25*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-0.5*pow(ma1,6) -
							0.5*pow(mpion,6) + pow(ma1,4)*(0.5*pow(mpion,2) +
							0.5*pow(mrho,2)) +
							               pow(mpion,4)*(0.5*pow(mrho,2) - 1.*s) +
							pow(mpion,2)*(-0.5*pow(mrho,2) + 0.5*s)*s +
							               pow(ma1,2)*(0.5*pow(mpion,4) +
							pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + (-0.5*pow(mrho,2) +
							0.5*s)*s)) +
							            eta1*(1.*pow(mpion,6) + pow(ma1,4)*(1.*pow(mpion,2) -
							0.5*s) + pow(mpion,2)*(1.5*pow(mrho,2) - 2.*s)*s +
							(-0.5*pow(mrho,2) + 0.5*s)*pow(s,2) +
							               pow(mpion,4)*(-1.*pow(mrho,2) + 1.5*s) +
							pow(ma1,2)*(-2.*pow(mpion,4) + 0.5*pow(mrho,2)*s +
							pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*
							          log(fabs(-1.*pow(ma1,2) + tmax)))/(1.*pow(ma1,2) -
							1.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s) -
							       (0.03125*pow(eta1 - 1.*eta2,2)*(1.*pow(ma1,2) -
							1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s)*
							          (eta1*eta2*(-2.*pow(ma1,8) +
							pow(ma1,6)*(8.*pow(mpion,2) + 4.*pow(mrho,2) - 4.*s) +
							               pow(ma1,4)*(-12.*pow(mpion,4) - 4.*pow(mrho,4) +
							4.*pow(mrho,2)*s + pow(mpion,2)*(-12.*pow(mrho,2) + 8.*s)) +
							               pow(mpion,2)*(-2.*pow(mpion,6) -
							4.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 4.*pow(mrho,4)*s
							- 2.*pow(mrho,2)*pow(s,2) +
							                  pow(mpion,2)*(-8.*pow(mrho,4) +
							8.*pow(mrho,2)*s)) +
							               pow(ma1,2)*(8.*pow(mpion,6) + 2.*pow(mrho,6) +
							pow(mpion,4)*(12.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,4)*s -
							2.*pow(mrho,2)*pow(s,2) + 2.*pow(s,3) +
							                  pow(mpion,2)*(8.*pow(mrho,4) - 4.*pow(mrho,2)*s
							- 4.*pow(s,2)))) +
							            pow(eta2,2)*(pow(ma1,8) +
							pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
							               pow(mpion,4)*(pow(mpion,4) +
							2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 1.*pow(mrho,2)*s) +
							               pow(ma1,4)*(6.*pow(mpion,4) - 1.*pow(mrho,4) +
							pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) + pow(mrho,2)*s) +
							               pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mrho,6) -
							5.*pow(mrho,4)*s + 4.*pow(mrho,2)*pow(s,2) - 1.*pow(s,3) +
							pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s) +
							                  pow(mpion,2)*(2.*pow(mrho,4) - 4.*pow(mrho,2)*s
							+ 2.*pow(s,2)))) +
							            pow(eta1,2)*(pow(ma1,8) + pow(mpion,8) +
							2.*pow(mpion,6)*pow(mrho,2) +
							pow(mpion,2)*pow(mrho,2)*s*(-2.*pow(mrho,2) + 2.*s) +
							               pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) +
							2.*s) + pow(mpion,4)*(3.*pow(mrho,4) - 5.*pow(mrho,2)*s) +
							               pow(ma1,4)*(6.*pow(mpion,4) + pow(mrho,4) +
							pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 3.*pow(mrho,2)*s) +
							               pow(ma1,2)*(-4.*pow(mpion,6) + pow(mrho,4)*s -
							1.*pow(s,3) + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s) +
							                  pow(mpion,2)*(-2.*pow(mrho,4) +
							4.*pow(mrho,2)*s + 2.*pow(s,2)))))*log(fabs(-1.*pow(ma1,2) + tmax)))/
							        (0.25*pow(Gammaa1,2)*pow(ma1,2) + 1.*pow(ma1,4) +
							1.*pow(mpion,4) + 1.*pow(mpion,2)*pow(mrho,2) + 0.25*pow(mrho,4)
							- 1.*pow(mpion,2)*s - 0.5*pow(mrho,2)*s +
							          0.25*pow(s,2) + pow(ma1,2)*(-2.*pow(mpion,2) -
							1.*pow(mrho,2) + 1.*s)) -
							       (1.*(1.*eta1 - 1.*eta2)*(eta2*(pow(ma1,4)*(0.5*pow(mrho,2) -
							1.*C4*pow(mrho,4)) +
							               pow(mpion,2)*pow(mrho,2)*(pow(mpion,2)*(0.5 -
							1.*C4*pow(mrho,2)) + (-0.25 + 0.125*delta)*(pow(mrho,2) + s)) +
							               pow(ma1,2)*(-1.*C4*pow(mrho,6) +
							pow(mpion,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4)) +
							0.25*delta*pow(s,2) + pow(mrho,2)*s*(-0.25 - 0.375*delta - 1.*C4*s) +
							                  pow(mrho,4)*(0.75 - 0.125*delta + 2.*C4*s))) +
							eta1*(pow(ma1,4)*(-0.5*pow(mrho,2) + 1.*C4*pow(mrho,4)) +
							               pow(ma1,2)*(0.5*pow(mrho,4) - 1.*C4*pow(mrho,6) +
							pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) -
							0.25*delta*pow(s,2) +
							                  1.*C4*pow(mrho,2)*pow(s,2)) +
							pow(mrho,2)*(pow(mpion,4)*(-0.5 + 1.*C4*pow(mrho,2)) + s*((0.25 -
							0.125*delta)*pow(mrho,2) + (-0.25 + 0.125*delta)*s) +
							                  pow(mpion,2)*(2.*C4*pow(mrho,4) + (0.5 +
							0.25*delta)*s + pow(mrho,2)*(-1. - 2.*C4*s)))))*log(fabs(-1.*pow(ma1,2) + tmax)))/pow(mrho,2) +
							       0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-1.*pow(mpion,2) + tmax)) +
							       (0.25*(0. + 8.000000000000002*pow(2. -
							1.*delta,2)*pow(mpion,4)*pow(mrho,2) - 5.999999999999999*pow(2. -
							1.*delta,2)*pow(mpion,2)*pow(mrho,2)*s +
							            1.*pow(2. -
							1.*delta,2)*pow(mrho,2)*pow(s,2))*log(fabs(-1.*pow(mpion,2) + tmax)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
							       (0.125*(-2. + delta)*(eta1 - 1.*eta2)*pow(mpion,2)*(0. +
							eta2*pow(mpion,2)*(4.*pow(mrho,2) - 4.*s) +
							            eta1*(2.*pow(mrho,4) - 2.*pow(mrho,2)*s +
							pow(mpion,2)*(-4.*pow(mrho,2) + 4.*s)))*log(fabs(-1.*pow(mpion,2) + tmax)))/(pow(ma1,2) - 1.*pow(mpion,2)) +
							       (2.*(-2. + 1.*delta)*(0. + (-0.25 + 0.125*delta)*pow(mrho,2)*s
							+ pow(mpion,2)*(-2.*C4*pow(mrho,4) - 0.5*delta*s + pow(mrho,2)*(1.
							+ 2.*C4*s)))*
							          log(fabs(-1.*pow(mpion,2) + tmax)))/pow(mrho,2) - (0.5*(-2. +
							delta)*(eta1 - 1.*eta2)*pow(mpion,2)*
							          (eta1*(pow(mpion,4)*(-1.*pow(mrho,2) + 1.*s) +
							pow(ma1,2)*(pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) +
							(-0.5*pow(mrho,2) + 0.5*s)*s) +
							               pow(mpion,2)*(-1.*pow(mrho,4) + 2.5*pow(mrho,2)*s
							- 1.5*pow(s,2)) + s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s +
							0.5*pow(s,2))) +
							            eta2*(0.5*pow(mrho,6) + pow(mpion,4)*(1.*pow(mrho,2)
							- 1.*s) - 1.5*pow(mrho,4)*s + 1.5*pow(mrho,2)*pow(s,2) -
							0.5*pow(s,3) +
							               pow(mpion,2)*(1.5*pow(mrho,4) - 3.*pow(mrho,2)*s +
							1.5*pow(s,2)) +
							               pow(ma1,2)*(-0.5*pow(mrho,4) + 1.*pow(mrho,2)*s -
							0.5*pow(s,2) + pow(mpion,2)*(-1.*pow(mrho,2) +
							1.*s))))*log(fabs(-1.*pow(mpion,2) + tmax)))/
							        (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4)
							+ 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 2.*pow(mpion,2)*s
							- 2.*pow(mrho,2)*s + pow(s,2) +
							          pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s))
							- 0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-1.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmax)) +
							       (0.5*(-2. + delta)*(eta1 -
							1.*eta2)*(eta2*pow(mpion,6)*(1.*pow(mrho,2) - 1.*s) +
							eta2*pow(ma1,2)*pow(mpion,4)*(-1.*pow(mrho,2) + 1.*s) +
							            eta1*pow(ma1,2)*pow(mpion,2)*(-0.5*pow(mrho,4) +
							pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + 0.5*pow(mrho,2)*s) +
							            eta1*pow(mpion,4)*(0.5*pow(mrho,4) -
							0.5*pow(mrho,2)*s + pow(mpion,2)*(-1.*pow(mrho,2) +
							1.*s)))*log(fabs(-1.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmax)))/
							        (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) -
							2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) +
							       (0.5*(-2. + delta)*(eta1 -
							1.*eta2)*pow(mpion,2)*(eta1*(pow(mpion,2)*(1.*pow(mrho,2) - 1.*s)
							+ (-0.5*pow(mrho,2) + 0.5*s)*s) +
							            eta2*(-0.5*pow(mrho,4) + 1.*pow(mrho,2)*s -
							0.5*pow(s,2) + pow(mpion,2)*(-1.*pow(mrho,2) +
							1.*s)))*log(fabs(-1.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmax)))/
							        (1.*pow(ma1,2) - 1.*pow(mpion,2) - 1.*pow(mrho,2) + 1.*s) -
							       (0.25*(0. + 8.000000000000002*pow(2. -
							1.*delta,2)*pow(mpion,4)*pow(mrho,2) + 1.*pow(2. -
							1.*delta,2)*pow(mrho,4)*s +
							            pow(mpion,2)*(C4*(32. - 16.*delta)*pow(mrho,6) +
							delta*(-8. + 4.*delta)*pow(s,2) +
							               pow(mrho,2)*s*(-8. + 24.*delta - 10.*pow(delta,2) +
							32.*C4*s - 16.*C4*delta*s) + pow(mrho,4)*(-16. + 8.*delta - 64.*C4*s +
							32.*C4*delta*s)))*
							          log(fabs(-1.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmax)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
							       0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(4.*pow(ma1,6) +
							pow(Gammaa1,2)*pow(ma1,2)*(-4.*pow(ma1,2) + 4.*pow(mpion,2) -
							2.*s) +
							             pow(ma1,4)*(-12.*pow(mpion,2) + 6.*s) +
							pow(mpion,2)*(-4.*pow(mpion,4) + 4.*pow(mrho,4) +
							2.*pow(mpion,2)*s - 2.*pow(mrho,2)*s) +
							             pow(ma1,2)*(12.*pow(mpion,4) - 2.*pow(mrho,4) -
							8.*pow(mpion,2)*s - 4.*pow(mrho,2)*s + 4.*pow(s,2))) +
							          pow(eta1,2)*(-2.*pow(ma1,6) + 2.*pow(mpion,6) +
							3.*pow(mpion,4)*pow(mrho,2) + pow(ma1,4)*(6.*pow(mpion,2) +
							3.*pow(mrho,2) - 3.*s) - 1.*pow(mpion,4)*s -
							             1.*pow(mpion,2)*pow(mrho,2)*s - 1.*pow(mrho,4)*s +
							pow(mrho,2)*pow(s,2) +
							             pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,2) -
							2.*pow(mpion,2) - 1.*pow(mrho,2) + s) +
							             pow(ma1,2)*(-6.*pow(mpion,4) - 1.*pow(mrho,4) +
							4.*pow(mrho,2)*s - 2.*pow(s,2) + pow(mpion,2)*(-6.*pow(mrho,2) +
							4.*s))) +
							          pow(eta2,2)*(-2.*pow(ma1,6) + 2.*pow(mpion,6) -
							3.*pow(mpion,4)*pow(mrho,2) + pow(ma1,4)*(6.*pow(mpion,2) -
							3.*pow(mrho,2) - 3.*s) - 1.*pow(mpion,4)*s +
							             pow(mpion,2)*pow(mrho,2)*s +
							pow(Gammaa1,2)*pow(ma1,2)*(2.*pow(ma1,2) - 2.*pow(mpion,2) +
							pow(mrho,2) + s) +
							             pow(ma1,2)*(-6.*pow(mpion,4) - 1.*pow(mrho,4) +
							2.*pow(mrho,2)*s - 2.*pow(s,2) + pow(mpion,2)*(6.*pow(mrho,2) +
							4.*s))))*
							        log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) -
							4.*pow(ma1,2)*pow(mpion,2) + 4.*pow(mpion,4) +
							2.*pow(ma1,2)*(-1.*pow(mrho,2) + s + tmax) -
							          4.*pow(mpion,2)*(-1.*pow(mrho,2) + s + tmax) +
							pow(-1.*pow(mrho,2) + s + tmax,2))) -
							       (0.5*(1.*eta1 -
							1.*eta2)*(eta2*(pow(Gammaa1,2)*pow(ma1,2)*pow(mrho,2)*(0.5 -
							1.*C4*pow(mrho,2)) + pow(ma1,4)*(-0.5*pow(mrho,2) +
							1.*C4*pow(mrho,4)) +
							               pow(mpion,2)*pow(mrho,2)*(pow(mpion,2)*(-0.5 +
							1.*C4*pow(mrho,2)) + (0.25 - 0.125*delta)*(pow(mrho,2) + s)) +
							               pow(ma1,2)*(1.*C4*pow(mrho,6) +
							pow(mpion,2)*(1.*pow(mrho,2) - 2.*C4*pow(mrho,4)) -
							0.25*delta*pow(s,2) + pow(mrho,4)*(-0.75 + 0.125*delta - 2.*C4*s) +
							                  pow(mrho,2)*s*(0.25 + 0.375*delta + 1.*C4*s))) + eta1*
							             (pow(Gammaa1,2)*pow(ma1,2)*pow(mrho,2)*(-0.5 +
							1.*C4*pow(mrho,2)) + pow(ma1,4)*(0.5*pow(mrho,2) -
							1.*C4*pow(mrho,4)) +
							               pow(ma1,2)*(-0.5*pow(mrho,4) + 1.*C4*pow(mrho,6) +
							pow(mpion,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4)) +
							0.25*delta*pow(s,2) -
							                  1.*C4*pow(mrho,2)*pow(s,2)) +
							pow(mrho,2)*(pow(mpion,4)*(0.5 - 1.*C4*pow(mrho,2)) + s*((-0.25 +
							0.125*delta)*pow(mrho,2) + (0.25 - 0.125*delta)*s) +
							                  pow(mpion,2)*(-2.*C4*pow(mrho,4) + (-0.5 -
							0.25*delta)*s + pow(mrho,2)*(1. + 2.*C4*s)))))*
							          log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(pow(ma1,2) -
							2.*pow(mpion,2) - 1.*pow(mrho,2) + s + tmax,2))))/pow(mrho,2) +
							       (0.125*(-2. + delta)*(eta1 -
							1.*eta2)*(eta2*(0.5*pow(Gammaa1,4)*pow(ma1,4) - 0.5*pow(ma1,8) +
							0.5*pow(mpion,8) + 0.5*pow(ma1,4)*pow(mpion,2)*pow(mrho,2) -
							               0.5*pow(mpion,6)*pow(mrho,2) +
							pow(Gammaa1,2)*(pow(ma1,2)*pow(mpion,2)*(1.*pow(mpion,2) +
							1.5*pow(mrho,2) - 2.*s) +
							                  pow(ma1,4)*(-1.*pow(mpion,2) + 0.5*pow(mrho,2)
							- 1.*s)) + pow(ma1,6)*(1.*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) +
							               pow(ma1,2)*pow(mpion,4)*(-1.*pow(mpion,2) -
							0.5*pow(mrho,2) + 1.*s)) +
							            eta1*(pow(ma1,6)*(1.*pow(mpion,2) + 0.5*s) +
							pow(ma1,2)*(3.*pow(mpion,6) + 1.*pow(mpion,2)*pow(mrho,4) -
							0.5*pow(mpion,4)*s) +
							               pow(ma1,4)*(-3.*pow(mpion,4) +
							pow(mpion,2)*(-1.*pow(mrho,2) + 0.5*s) - 0.5*pow(mrho,2)*s) +
							               pow(mpion,4)*(-1.*pow(mpion,4) - 1.*pow(mrho,4) +
							pow(mpion,2)*(1.*pow(mrho,2) - 0.5*s) + 0.5*pow(mrho,2)*s) +
							               pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(mpion,4) +
							pow(ma1,2)*(1.*pow(mpion,2) + 0.5*s) - 0.5*pow(mrho,2)*s +
							pow(mpion,2)*(-1.*pow(mrho,2) + 1.5*s))))*
							          log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) +
							4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
							4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s +
							            pow(s,2) - 4.*pow(mpion,2)*tmax - 2.*pow(mrho,2)*tmax +
							2.*s*tmax + pow(tmax,2) + pow(ma1,2)*(-4.*pow(mpion,2) -
							2.*pow(mrho,2) + 2.*s + 2.*tmax))))/
							        (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) -
							2.*pow(ma1,2)*pow(mpion,2) + pow(mpion,4)) -
							       (0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta1*(-1.*pow(mpion,8)
							+ 1.*pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(1.*pow(mpion,2) -
							0.5*s) - 0.5*pow(mpion,6)*s -
							               4.*pow(mpion,4)*pow(mrho,2)*s -
							1.5*pow(mpion,2)*pow(mrho,4)*s + 3.5*pow(mpion,4)*pow(s,2) +
							4.*pow(mpion,2)*pow(mrho,2)*pow(s,2) +
							               0.5*pow(mrho,4)*pow(s,2) -
							2.5*pow(mpion,2)*pow(s,3) - 1.*pow(mrho,2)*pow(s,3) +
							0.5*pow(s,4) +
							               pow(Gammaa1,2)*pow(ma1,2)*(-1.*pow(mpion,4) +
							pow(ma1,2)*(1.*pow(mpion,2) - 0.5*s) - 0.5*pow(mpion,2)*s +
							0.5*pow(s,2)) +
							               pow(ma1,4)*(-3.*pow(mpion,4) + (1.*pow(mrho,2) -
							0.5*s)*s + pow(mpion,2)*(-2.*pow(mrho,2) + 2.5*s)) +
							               pow(ma1,2)*(3.*pow(mpion,6) +
							pow(mpion,4)*(2.*pow(mrho,2) - 1.5*s) - 0.5*pow(mrho,4)*s +
							0.5*pow(s,3) +
							                  pow(mpion,2)*(1.*pow(mrho,4) - 1.*pow(mrho,2)*s
							- 1.*pow(s,2)))) +
							            eta2*(0.5*pow(Gammaa1,4)*pow(ma1,4) - 0.5*pow(ma1,8)
							+ 0.5*pow(mpion,8) + pow(ma1,6)*(1.*pow(mpion,2) +
							1.*pow(mrho,2) - 0.5*s) + 0.5*pow(mpion,6)*s +
							               pow(ma1,4)*(-0.5*pow(mrho,4) + (-0.5*pow(mpion,2)
							+ 0.5*s)*s) + pow(mpion,4)*(-0.5*pow(mrho,4) + 2.*pow(mrho,2)*s -
							1.5*pow(s,2)) +
							               pow(mpion,2)*s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s
							+ 0.5*pow(s,2)) +
							               pow(Gammaa1,2)*pow(ma1,2)*(1.*pow(mpion,4) +
							0.5*pow(mrho,4) + pow(mpion,2)*(2.*pow(mrho,2) - 1.5*s) -
							1.*pow(mrho,2)*s + 0.5*pow(s,2) +
							                  pow(ma1,2)*(-1.*pow(mpion,2) - 1.*pow(mrho,2) +
							1.5*s)) +
							               pow(ma1,2)*(-1.*pow(mpion,6) +
							pow(mpion,4)*(-1.*pow(mrho,2) + 0.5*s) +
							pow(mpion,2)*(-1.*pow(mrho,4) + 2.*pow(mrho,2)*s - 1.*pow(s,2)) +
							                  s*(0.5*pow(mrho,4) - 1.*pow(mrho,2)*s +
							0.5*pow(s,2)))))*
							          log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) +
							4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
							4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s +
							            pow(s,2) - 4.*pow(mpion,2)*tmax - 2.*pow(mrho,2)*tmax +
							2.*s*tmax + pow(tmax,2) + pow(ma1,2)*(-4.*pow(mpion,2) -
							2.*pow(mrho,2) + 2.*s + 2.*tmax))))/
							        (pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) + pow(mpion,4)
							+ 2.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) - 2.*pow(mpion,2)*s
							- 2.*pow(mrho,2)*s + pow(s,2) +
							          pow(ma1,2)*(-2.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s))
							- (0.0625*pow(eta1 - 1.*eta2,2)*
							          (pow(eta2,2)*(-1.*pow(ma1,10) +
							pow(ma1,8)*(5.*pow(mpion,2) + 2.5*pow(mrho,2) - 2.5*s) +
							               pow(Gammaa1,4)*pow(ma1,4)*(1.*pow(ma1,2) -
							1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s) +
							               pow(ma1,4)*(10.*pow(mpion,6) - 2.5*pow(mrho,6) +
							pow(mpion,4)*(15.*pow(mrho,2) - 9.*s) + 6.*pow(mrho,4)*s -
							4.5*pow(mrho,2)*pow(s,2) + 1.*pow(s,3)) +
							               pow(ma1,6)*(-10.*pow(mpion,4) + (1.*pow(mrho,2) -
							1.*s)*s + pow(mpion,2)*(-10.*pow(mrho,2) + 8.*s)) +
							               pow(mpion,4)*(1.*pow(mpion,6) + 0.5*pow(mrho,6) +
							pow(mpion,4)*(2.5*pow(mrho,2) - 0.5*s) - 1.*pow(mrho,4)*s +
							0.5*pow(mrho,2)*pow(s,2) +
							                  pow(mpion,2)*(2.*pow(mrho,4) -
							2.*pow(mrho,2)*s)) +
							               pow(Gammaa1,2)*pow(ma1,2)*(4.*pow(ma1,6) -
							4.*pow(mpion,6) - 0.5*pow(mrho,6) + 1.5*pow(mrho,4)*s -
							1.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) +
							                  pow(mpion,4)*(-6.*pow(mrho,2) + 6.*s) +
							pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
							                  pow(mpion,2)*(-3.*pow(mrho,4) +
							6.*pow(mrho,2)*s - 3.*pow(s,2)) +
							                  pow(ma1,2)*(12.*pow(mpion,4) + 3.*pow(mrho,4) +
							pow(mpion,2)*(12.*pow(mrho,2) - 12.*s) - 6.*pow(mrho,2)*s +
							3.*pow(s,2))) +
							               pow(ma1,2)*(-5.*pow(mpion,8) + 1.*pow(mrho,8) -
							3.5*pow(mrho,6)*s + 4.5*pow(mrho,4)*pow(s,2) -
							2.5*pow(mrho,2)*pow(s,3) + 0.5*pow(s,4) +
							                  pow(mpion,6)*(-10.*pow(mrho,2) + 4.*s) +
							pow(mpion,4)*(-2.*pow(mrho,4) + 1.*pow(mrho,2)*s + 1.*pow(s,2)) +
							                  pow(mpion,2)*(3.*pow(mrho,6) - 8.*pow(mrho,4)*s
							+ 7.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3)))) +
							            pow(eta1,2)*(-1.*pow(ma1,10) +
							pow(ma1,8)*(5.*pow(mpion,2) + 2.5*pow(mrho,2) - 2.5*s) +
							               pow(Gammaa1,4)*pow(ma1,4)*(1.*pow(ma1,2) -
							1.*pow(mpion,2) - 0.5*pow(mrho,2) + 0.5*s) +
							               pow(ma1,6)*(-10.*pow(mpion,4) - 2.*pow(mrho,4) +
							5.*pow(mrho,2)*s - 1.*pow(s,2) + pow(mpion,2)*(-10.*pow(mrho,2)
							+ 8.*s)) +
							               pow(ma1,4)*(10.*pow(mpion,6) + 0.5*pow(mrho,6) +
							pow(mpion,4)*(15.*pow(mrho,2) - 9.*s) - 3.*pow(mrho,4)*s +
							1.5*pow(mrho,2)*pow(s,2) + 1.*pow(s,3) +
							                  pow(mpion,2)*(6.*pow(mrho,4) -
							12.*pow(mrho,2)*s)) +
							               pow(Gammaa1,2)*pow(ma1,2)*(4.*pow(ma1,6) -
							4.*pow(mpion,6) - 0.5*pow(mrho,6) + 1.5*pow(mrho,4)*s -
							1.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3) +
							                  pow(mpion,4)*(-6.*pow(mrho,2) + 6.*s) +
							pow(ma1,4)*(-12.*pow(mpion,2) - 6.*pow(mrho,2) + 6.*s) +
							                  pow(mpion,2)*(-3.*pow(mrho,4) +
							6.*pow(mrho,2)*s - 3.*pow(s,2)) +
							                  pow(ma1,2)*(12.*pow(mpion,4) + 3.*pow(mrho,4) +
							pow(mpion,2)*(12.*pow(mrho,2) - 12.*s) - 6.*pow(mrho,2)*s +
							3.*pow(s,2))) +
							               pow(mpion,2)*(1.*pow(mpion,8) +
							pow(mpion,6)*(2.5*pow(mrho,2) - 0.5*s) +
							pow(mpion,4)*(4.*pow(mrho,4) - 6.*pow(mrho,2)*s) +
							                  pow(mrho,2)*s*(-1.*pow(mrho,4) +
							2.*pow(mrho,2)*s - 1.*pow(s,2)) + pow(mpion,2)*(1.5*pow(mrho,6)
							- 6.*pow(mrho,4)*s + 4.5*pow(mrho,2)*pow(s,2))) +
							               pow(ma1,2)*(-5.*pow(mpion,8) +
							pow(mpion,6)*(-10.*pow(mrho,2) + 4.*s) +
							pow(mpion,4)*(-8.*pow(mrho,4) + 13.*pow(mrho,2)*s + 1.*pow(s,2)) +
							                  pow(mpion,2)*(-1.*pow(mrho,6) +
							6.*pow(mrho,4)*s - 3.*pow(mrho,2)*pow(s,2) - 2.*pow(s,3)) +
							                  s*(0.5*pow(mrho,6) - 0.5*pow(mrho,4)*s -
							0.5*pow(mrho,2)*pow(s,2) + 0.5*pow(s,3)))) +
							            eta1*eta2*(2.*pow(ma1,10) +
							pow(Gammaa1,4)*pow(ma1,4)*(-2.*pow(ma1,2) + 2.*pow(mpion,2) +
							1.*pow(mrho,2) - 1.*s) +
							               pow(ma1,8)*(-10.*pow(mpion,2) - 5.*pow(mrho,2) +
							5.*s) +
							               pow(ma1,6)*(20.*pow(mpion,4) + 6.*pow(mrho,4) +
							pow(mpion,2)*(20.*pow(mrho,2) - 16.*s) - 8.*pow(mrho,2)*s +
							2.*pow(s,2)) +
							               pow(ma1,4)*(-20.*pow(mpion,6) - 4.*pow(mrho,6) +
							6.*pow(mrho,4)*s - 2.*pow(s,3) + pow(mpion,4)*(-30.*pow(mrho,2)
							+ 18.*s) +
							                  pow(mpion,2)*(-18.*pow(mrho,4) +
							18.*pow(mrho,2)*s)) +
							               pow(mpion,2)*(-2.*pow(mpion,8) - 1.*pow(mrho,8) +
							3.*pow(mrho,6)*s - 3.*pow(mrho,4)*pow(s,2) +
							1.*pow(mrho,2)*pow(s,3) +
							                  pow(mpion,6)*(-5.*pow(mrho,2) + 1.*s) +
							pow(mpion,4)*(-10.*pow(mrho,4) + 10.*pow(mrho,2)*s) +
							                  pow(mpion,2)*(-6.*pow(mrho,6) +
							12.*pow(mrho,4)*s - 6.*pow(mrho,2)*pow(s,2))) +
							               pow(ma1,2)*(10.*pow(mpion,8) + 1.*pow(mrho,8) +
							pow(mpion,6)*(20.*pow(mrho,2) - 8.*s) - 2.*pow(mrho,6)*s +
							2.*pow(mrho,2)*pow(s,3) - 1.*pow(s,4) +
							                  pow(mpion,4)*(22.*pow(mrho,4) -
							20.*pow(mrho,2)*s - 2.*pow(s,2)) + pow(mpion,2)*(8.*pow(mrho,6)
							- 12.*pow(mrho,4)*s + 4.*pow(s,3))) +
							               pow(Gammaa1,2)*pow(ma1,2)*(-8.*pow(ma1,6) +
							8.*pow(mpion,6) + 1.*pow(mrho,6) + pow(mpion,4)*(12.*pow(mrho,2)
							- 12.*s) +
							                  pow(ma1,4)*(24.*pow(mpion,2) + 12.*pow(mrho,2)
							- 12.*s) - 3.*pow(mrho,4)*s + 3.*pow(mrho,2)*pow(s,2) -
							1.*pow(s,3) +
							                  pow(mpion,2)*(6.*pow(mrho,4) -
							12.*pow(mrho,2)*s + 6.*pow(s,2)) +
							                  pow(ma1,2)*(-24.*pow(mpion,4) - 6.*pow(mrho,4)
							+ 12.*pow(mrho,2)*s - 6.*pow(s,2) +
							pow(mpion,2)*(-24.*pow(mrho,2) + 24.*s)))))*
							          log(fabs(pow(Gammaa1,2)*pow(ma1,2) + pow(ma1,4) +
							4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
							4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s +
							            pow(s,2) - 4.*pow(mpion,2)*tmax - 2.*pow(mrho,2)*tmax +
							2.*s*tmax + pow(tmax,2) + pow(ma1,2)*(-4.*pow(mpion,2) -
							2.*pow(mrho,2) + 2.*s + 2.*tmax))))/
							        (pow(Gammaa1,2)*pow(ma1,2) + 4.*pow(ma1,4) +
							4.*pow(mpion,4) + 4.*pow(mpion,2)*pow(mrho,2) + pow(mrho,4) -
							4.*pow(mpion,2)*s - 2.*pow(mrho,2)*s + pow(s,2) +
							          pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) +
							4.*s))))/(16.*Pi*s*(-4*pow(mpion,2) + s)));

  return sigma / spin_deg_factor;
}

//C22
double total_xsection_C22(double s) {
  double sigma;
  double spin_deg_factor = 1.0;
  double tmin = min_mandelstam_t(s, mpion, mpion, mrho);
  double tmax = max_mandelstam_t(s, mpion, mpion, mrho);

  sigma = (pow(Const,2)*pow(ghat,4)*((-0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*
	             (-2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) +
	               pow(ma1,2)*pow(mpion,2)*(8.*pow(mpion,4) - 8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
	               pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
	            pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) + 1.*pow(mpion,4)*pow(mrho,4) +
	               pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
	               pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2)) +
	               pow(ma1,2)*(-4.*pow(mpion,6) - 2.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(6.*pow(mrho,2) + 2.*s))) +
	            pow(eta1,2)*(1.*pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
	               pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2)) +
	               pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
	               pow(mpion,2)*(1.*pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s + pow(mpion,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s)))))/
	        (1.*pow(ma1,2) - 1.*tmin) - (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/(1.*pow(mpion,2) - 1.*tmin) - 0.75*tmin +
	       (0.25*pow(-2. + delta,2)*pow(mpion,2)*tmin)/pow(mrho,2) + 0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-1.*pow(ma1,2) + pow(mrho,2) - 2.*s) + eta1*(2.*pow(mpion,2) + s))*
	        tmin - C4*(2.*C4*pow(mrho,4) + pow(mrho,2)*(-3. - 4.*C4*s) + s*(3. + 2.*C4*s))*tmin -
	       (0.5*pow(delta,2)*(1.*pow(mpion,4)*pow(mrho,2) + 0.25*pow(mrho,6) - 0.75*pow(mrho,4)*s + 0.125*pow(mrho,2)*pow(s,2) + 0.25*pow(s,3) +
	            pow(mpion,2)*(2.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.75*pow(s,2)))*tmin)/pow(mrho,6) -
	       0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-6.*pow(ma1,4) - 12.*pow(mpion,4) + 2.*pow(mrho,4) + pow(ma1,2)*(16.*pow(mpion,2) - 8.*s) + 8.*pow(mpion,2)*s +
	             4.*pow(mrho,2)*s - 4.*pow(s,2)) + pow(eta1,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
	             2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
	          pow(eta2,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2) +
	             pow(ma1,2)*(-8.*pow(mpion,2) + 4.*pow(mrho,2) + 4.*s)))*tmin -
	       (3.*delta*(0.6666666666666666*C4*pow(mrho,6) - 0.08333333333333333*pow(s,2) + pow(mrho,4)*(-0.25 - 0.5*C4*s) +
	            pow(mrho,2)*s*(0.08333333333333333 - 0.16666666666666666*C4*s) +
	            pow(mpion,2)*(1.*C4*pow(mrho,4) + 0.08333333333333333*s + pow(mrho,2)*(-0.4166666666666667 - 0.3333333333333333*C4*s)))*tmin)/pow(mrho,4) -
	       (0.25*(1.*eta1 - 1.*eta2)*(pow(mrho,2)*(eta1*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(ma1,2)*(1. - 2.*C4*pow(mrho,2)) + pow(mpion,2)*(-2. + 4.*C4*pow(mrho,2)) -
	                  2.*C4*pow(s,2)) + eta2*(-1.5*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) + pow(ma1,2)*(-1. + 2.*C4*pow(mrho,2)) + 0.5*s -
	                  4.*C4*pow(mrho,2)*s + 2.*C4*pow(s,2))) + delta*(eta2*
	                (-1.*pow(ma1,4) - 3.*pow(mpion,4) + 1.*pow(mrho,4) + pow(ma1,2)*(3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) + 0.25*pow(mrho,2)*s - 0.75*pow(s,2) +
	                  pow(mpion,2)*(1.*pow(mrho,2) + 1.*s)) + eta1*(1.*pow(ma1,4) + 3.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) -
	                  2.*pow(mrho,2)*s + 1.*pow(s,2) + pow(ma1,2)*(-3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s))))*tmin)/pow(mrho,2) -
	       0.0625*(-2. + delta)*(eta1 - 1.*eta2)*eta2*pow(tmin,2) + (0.25*pow(delta,2)*(2.*pow(mpion,2)*pow(mrho,2) + 1.*pow(mrho,4) - 0.75*pow(mrho,2)*s - 0.25*pow(s,2))*
	          pow(tmin,2))/pow(mrho,6) - 0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) + 2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
	          eta1*(pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*pow(tmin,2) +
	       (1.5*delta*(1.*C4*pow(mrho,4) + 0.08333333333333333*s + pow(mrho,2)*(-0.4166666666666667 - 0.3333333333333333*C4*s))*pow(tmin,2))/pow(mrho,4) -
	       (0.125*(1.*eta1 - 1.*eta2)*(pow(mrho,2)*(eta1*(1. - 2.*C4*pow(mrho,2)) + eta2*(-1. + 2.*C4*pow(mrho,2))) +
	            delta*(eta2*(-1.*pow(ma1,2) + 3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) + eta1*(1.*pow(ma1,2) - 3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s)))*pow(tmin,2))/
	        pow(mrho,2) - 0.010416666666666666*pow(eta1 - 1.*eta2,4)*pow(tmin,3) - (0.16666666666666666*pow(delta,2)*pow(tmin,3))/pow(mrho,4) -
	       (0.08333333333333333*delta*pow(1.*eta1 - 1.*eta2,2)*pow(tmin,3))/pow(mrho,2) +
	       (0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-2.*pow(ma1,8) - 2.*pow(mpion,8) + 2.*pow(mpion,4)*pow(mrho,4) + pow(ma1,6)*(8.*pow(mpion,2) - 4.*s) +
	               pow(ma1,2)*pow(mpion,2)*(8.*pow(mpion,4) - 8.*pow(mrho,4) - 4.*pow(mpion,2)*s + 4.*pow(mrho,2)*s) +
	               pow(ma1,4)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
	            pow(eta2,2)*(1.*pow(ma1,8) + 1.*pow(mpion,8) - 2.*pow(mpion,6)*pow(mrho,2) + 1.*pow(mpion,4)*pow(mrho,4) +
	               pow(ma1,6)*(-4.*pow(mpion,2) + 2.*pow(mrho,2) + 2.*s) +
	               pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2)) +
	               pow(ma1,2)*(-4.*pow(mpion,6) - 2.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(6.*pow(mrho,2) + 2.*s))) +
	            pow(eta1,2)*(1.*pow(ma1,8) + pow(ma1,6)*(-4.*pow(mpion,2) - 2.*pow(mrho,2) + 2.*s) +
	               pow(ma1,4)*(6.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2)) +
	               pow(ma1,2)*(-4.*pow(mpion,6) + 2.*pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(2.*pow(mrho,2) - 2.*s)*s + pow(mpion,4)*(-6.*pow(mrho,2) + 2.*s)) +
	               pow(mpion,2)*(1.*pow(mpion,6) + 2.*pow(mpion,4)*pow(mrho,2) - 2.*pow(mrho,6) + 2.*pow(mrho,4)*s + pow(mpion,2)*(1.*pow(mrho,4) - 2.*pow(mrho,2)*s)))))/
	        (1.*pow(ma1,2) - 1.*tmax) + (1.*pow(-2. + delta,2)*pow(mpion,2)*(1.*pow(mpion,2) - 0.25*pow(mrho,2)))/(1.*pow(mpion,2) - 1.*tmax) + 0.75*tmax -
	       (0.25*pow(-2. + delta,2)*pow(mpion,2)*tmax)/pow(mrho,2) - 0.125*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-1.*pow(ma1,2) + pow(mrho,2) - 2.*s) + eta1*(2.*pow(mpion,2) + s))*
	        tmax + C4*(2.*C4*pow(mrho,4) + pow(mrho,2)*(-3. - 4.*C4*s) + s*(3. + 2.*C4*s))*tmax +
	       (0.5*pow(delta,2)*(1.*pow(mpion,4)*pow(mrho,2) + 0.25*pow(mrho,6) - 0.75*pow(mrho,4)*s + 0.125*pow(mrho,2)*pow(s,2) + 0.25*pow(s,3) +
	            pow(mpion,2)*(2.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.75*pow(s,2)))*tmax)/pow(mrho,6) +
	       0.03125*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-6.*pow(ma1,4) - 12.*pow(mpion,4) + 2.*pow(mrho,4) + pow(ma1,2)*(16.*pow(mpion,2) - 8.*s) + 8.*pow(mpion,2)*s +
	             4.*pow(mrho,2)*s - 4.*pow(s,2)) + pow(eta1,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s +
	             2.*pow(s,2) + pow(ma1,2)*(-8.*pow(mpion,2) - 4.*pow(mrho,2) + 4.*s)) +
	          pow(eta2,2)*(3.*pow(ma1,4) + 6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2) +
	             pow(ma1,2)*(-8.*pow(mpion,2) + 4.*pow(mrho,2) + 4.*s)))*tmax +
	       (3.*delta*(0.6666666666666666*C4*pow(mrho,6) - 0.08333333333333333*pow(s,2) + pow(mrho,4)*(-0.25 - 0.5*C4*s) +
	            pow(mrho,2)*s*(0.08333333333333333 - 0.16666666666666666*C4*s) +
	            pow(mpion,2)*(1.*C4*pow(mrho,4) + 0.08333333333333333*s + pow(mrho,2)*(-0.4166666666666667 - 0.3333333333333333*C4*s)))*tmax)/pow(mrho,4) +
	       (0.25*(1.*eta1 - 1.*eta2)*(pow(mrho,2)*(eta1*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(ma1,2)*(1. - 2.*C4*pow(mrho,2)) + pow(mpion,2)*(-2. + 4.*C4*pow(mrho,2)) -
	                  2.*C4*pow(s,2)) + eta2*(-1.5*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) + pow(ma1,2)*(-1. + 2.*C4*pow(mrho,2)) + 0.5*s -
	                  4.*C4*pow(mrho,2)*s + 2.*C4*pow(s,2))) + delta*(eta2*
	                (-1.*pow(ma1,4) - 3.*pow(mpion,4) + 1.*pow(mrho,4) + pow(ma1,2)*(3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) + 0.25*pow(mrho,2)*s - 0.75*pow(s,2) +
	                  pow(mpion,2)*(1.*pow(mrho,2) + 1.*s)) + eta1*(1.*pow(ma1,4) + 3.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) -
	                  2.*pow(mrho,2)*s + 1.*pow(s,2) + pow(ma1,2)*(-3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s))))*tmax)/pow(mrho,2) +
	       0.0625*(-2. + delta)*(eta1 - 1.*eta2)*eta2*pow(tmax,2) - (0.25*pow(delta,2)*(2.*pow(mpion,2)*pow(mrho,2) + 1.*pow(mrho,4) - 0.75*pow(mrho,2)*s - 0.25*pow(s,2))*
	          pow(tmax,2))/pow(mrho,6) + 0.03125*pow(eta1 - 1.*eta2,3)*(eta2*(-1.*pow(ma1,2) + 2.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
	          eta1*(pow(ma1,2) - 2.*pow(mpion,2) - 1.*pow(mrho,2) + s))*pow(tmax,2) -
	       (1.5*delta*(1.*C4*pow(mrho,4) + 0.08333333333333333*s + pow(mrho,2)*(-0.4166666666666667 - 0.3333333333333333*C4*s))*pow(tmax,2))/pow(mrho,4) +
	       (0.125*(1.*eta1 - 1.*eta2)*(pow(mrho,2)*(eta1*(1. - 2.*C4*pow(mrho,2)) + eta2*(-1. + 2.*C4*pow(mrho,2))) +
	            delta*(eta2*(-1.*pow(ma1,2) + 3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) + eta1*(1.*pow(ma1,2) - 3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s)))*pow(tmax,2))/
	        pow(mrho,2) + 0.010416666666666666*pow(eta1 - 1.*eta2,4)*pow(tmax,3) + (0.16666666666666666*pow(delta,2)*pow(tmax,3))/pow(mrho,4) +
	       (0.08333333333333333*delta*pow(1.*eta1 - 1.*eta2,2)*pow(tmax,3))/pow(mrho,2) - (2.*pow(mpion,4)*tmin*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) -
	       (5.*pow(mpion,2)*pow(mrho,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) -
	       (0.5*pow(mrho,4)*tmin*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) - (2.499999*pow(mpion,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (5.0000100000000005*delta*pow(mpion,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (1.2500010000000001*pow(mrho,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (0.500001*delta*pow(mrho,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (6.*C4*pow(mpion,2)*pow(mrho,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) + (3.*C4*pow(mrho,4)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (5.*delta*pow(mpion,2)*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) +
	       (1.75*pow(mrho,2)*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) +
	       (0.5*delta*pow(mrho,2)*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) + (1.5*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (1.0000005*delta*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) - (0.2500005*pow(delta,2)*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (2.000001*C4*pow(mpion,2)*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (3.*C4*delta*pow(mpion,2)*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (3.9999899999999995*C4*pow(mrho,2)*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (1.5*C4*delta*pow(mrho,2)*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (1.75*delta*pow(s,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) -
	       (0.125*pow(delta,2)*pow(s,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) -
	       (0.5*pow(delta,2)*pow(mpion,4)*pow(s,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) +
	       (0.999999*C4*pow(s,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (1.9999949999999997*C4*delta*pow(s,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (1.*delta*pow(mpion,2)*pow(s,3)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) +
	       (0.25*pow(delta,2)*pow(mpion,2)*pow(s,4)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6)*pow(pow(mrho,2) - 1.*s,2)) +
	       (0.25*delta*pow(s,4)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) -
	       (0.0625*pow(delta,2)*pow(s,5)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6)*pow(pow(mrho,2) - 1.*s,2)) +
	       (2.*delta*pow(mpion,4)*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) +
	       (1.*pow(mpion,2)*pow(s,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) -
	       (1.25*pow(delta,2)*pow(mpion,2)*pow(s,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) -
	       (0.25*pow(s,3)*tmin*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) +
	       (0.4375*pow(delta,2)*pow(s,3)*tmin*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) +
	       (2.000001*delta*pow(mpion,4)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (0.500001*pow(mpion,2)*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (1.4999993999999999*delta*pow(mpion,2)*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (2.5000050000000003*pow(delta,2)*pow(mpion,2)*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (0.2499999*pow(s,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (0.9999998999999999*delta*pow(s,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (0.8125005000000001*pow(delta,2)*pow(s,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (1.0000005*C4*delta*pow(mpion,2)*pow(s,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (0.4999995*C4*delta*pow(s,3)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (1.0000005*pow(delta,2)*pow(mpion,4)*s*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) -
	       (1.5000015000000002*delta*pow(mpion,2)*pow(s,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) -
	       (0.12499995*pow(delta,2)*pow(mpion,2)*pow(s,2)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) +
	       (0.49999994999999997*delta*pow(s,3)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) +
	       (0.12499995*pow(delta,2)*pow(s,3)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) +
	       (0.6250005000000001*pow(delta,2)*pow(mpion,2)*pow(s,3)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,8) - 1.*pow(mrho,6)*s) -
	       (0.1875*pow(delta,2)*pow(s,4)*tmin*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,8) - 1.*pow(mrho,6)*s) +
	       (2.*pow(mpion,2)*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) +
	       (1.*pow(mrho,2)*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) + (1.2499995*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (1.0000005*delta*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) - (3.*C4*pow(mrho,2)*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (1.*s*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) - (1.*delta*s*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) +
	       (1.0000005*C4*s*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) + (1.5*C4*delta*s*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (0.5*pow(delta,2)*pow(mpion,2)*pow(s,2)*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) -
	       (0.25*pow(delta,2)*pow(s,3)*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) -
	       (2.*delta*pow(mpion,2)*s*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) +
	       (1.*delta*pow(s,2)*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) +
	       (0.25*pow(delta,2)*pow(s,2)*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) -
	       (1.9999949999999997*delta*pow(mpion,2)*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (0.2500005*s*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (0.24999974999999997*delta*s*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (0.50000025*pow(delta,2)*s*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (0.50000025*C4*delta*pow(s,2)*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (0.9999974999999999*pow(delta,2)*pow(mpion,2)*s*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) +
	       (0.2500002*delta*pow(s,2)*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) -
	       (0.43749974999999997*pow(delta,2)*pow(s,2)*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) -
	       (0.062499975*pow(delta,2)*pow(s,3)*pow(tmin,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,8) - 1.*pow(mrho,6)*s) -
	       (0.6666666666666666*pow(tmin,3)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) -
	       (0.16666666666666666*pow(delta,2)*pow(s,2)*pow(tmin,3)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) +
	       (0.6666666666666666*delta*s*pow(tmin,3)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) +
	       (0.666667*delta*pow(tmin,3)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (0.3333335*pow(delta,2)*s*pow(tmin,3)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) +
	       (2.*pow(mpion,4)*tmax*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) +
	       (5.*pow(mpion,2)*pow(mrho,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) +
	       (0.5*pow(mrho,4)*tmax*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) + (2.499999*pow(mpion,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (5.0000100000000005*delta*pow(mpion,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (1.2500010000000001*pow(mrho,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (0.500001*delta*pow(mrho,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (6.*C4*pow(mpion,2)*pow(mrho,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) - (3.*C4*pow(mrho,4)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (5.*delta*pow(mpion,2)*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) -
	       (1.75*pow(mrho,2)*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) -
	       (0.5*delta*pow(mrho,2)*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) - (1.5*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (1.0000005*delta*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) + (0.2500005*pow(delta,2)*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (2.000001*C4*pow(mpion,2)*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (3.*C4*delta*pow(mpion,2)*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (3.9999899999999995*C4*pow(mrho,2)*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (1.5*C4*delta*pow(mrho,2)*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (1.75*delta*pow(s,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) +
	       (0.125*pow(delta,2)*pow(s,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) +
	       (0.5*pow(delta,2)*pow(mpion,4)*pow(s,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) -
	       (0.999999*C4*pow(s,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (1.9999949999999997*C4*delta*pow(s,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (1.*delta*pow(mpion,2)*pow(s,3)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) -
	       (0.25*pow(delta,2)*pow(mpion,2)*pow(s,4)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6)*pow(pow(mrho,2) - 1.*s,2)) -
	       (0.25*delta*pow(s,4)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) +
	       (0.0625*pow(delta,2)*pow(s,5)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6)*pow(pow(mrho,2) - 1.*s,2)) -
	       (2.*delta*pow(mpion,4)*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) -
	       (1.*pow(mpion,2)*pow(s,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) +
	       (1.25*pow(delta,2)*pow(mpion,2)*pow(s,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) +
	       (0.25*pow(s,3)*tmax*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) -
	       (0.4375*pow(delta,2)*pow(s,3)*tmax*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) -
	       (2.000001*delta*pow(mpion,4)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (0.500001*pow(mpion,2)*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (1.4999993999999999*delta*pow(mpion,2)*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (2.5000050000000003*pow(delta,2)*pow(mpion,2)*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (0.2499999*pow(s,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (0.9999998999999999*delta*pow(s,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (0.8125005000000001*pow(delta,2)*pow(s,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (1.0000005*C4*delta*pow(mpion,2)*pow(s,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (0.4999995*C4*delta*pow(s,3)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (1.0000005*pow(delta,2)*pow(mpion,4)*s*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) +
	       (1.5000015000000002*delta*pow(mpion,2)*pow(s,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) +
	       (0.12499995*pow(delta,2)*pow(mpion,2)*pow(s,2)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) -
	       (0.49999994999999997*delta*pow(s,3)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) -
	       (0.12499995*pow(delta,2)*pow(s,3)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) -
	       (0.6250005000000001*pow(delta,2)*pow(mpion,2)*pow(s,3)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,8) - 1.*pow(mrho,6)*s) +
	       (0.1875*pow(delta,2)*pow(s,4)*tmax*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,8) - 1.*pow(mrho,6)*s) -
	       (2.*pow(mpion,2)*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) -
	       (1.*pow(mrho,2)*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) - (1.2499995*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (1.0000005*delta*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) + (3.*C4*pow(mrho,2)*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) +
	       (1.*s*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) + (1.*delta*s*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) -
	       (1.0000005*C4*s*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) - (1.5*C4*delta*s*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,2) - 1.*s) -
	       (0.5*pow(delta,2)*pow(mpion,2)*pow(s,2)*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) +
	       (0.25*pow(delta,2)*pow(s,3)*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) +
	       (2.*delta*pow(mpion,2)*s*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) -
	       (1.*delta*pow(s,2)*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) -
	       (0.25*pow(delta,2)*pow(s,2)*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) +
	       (1.9999949999999997*delta*pow(mpion,2)*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (0.2500005*s*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (0.24999974999999997*delta*s*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (0.50000025*pow(delta,2)*s*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (0.50000025*C4*delta*pow(s,2)*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) -
	       (0.9999974999999999*pow(delta,2)*pow(mpion,2)*s*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) -
	       (0.2500002*delta*pow(s,2)*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) +
	       (0.43749974999999997*pow(delta,2)*pow(s,2)*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) +
	       (0.062499975*pow(delta,2)*pow(s,3)*pow(tmax,2)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,8) - 1.*pow(mrho,6)*s) +
	       (0.6666666666666666*pow(tmax,3)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,2) - 1.*s,2) +
	       (0.16666666666666666*pow(delta,2)*pow(s,2)*pow(tmax,3)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4)*pow(pow(mrho,2) - 1.*s,2)) -
	       (0.6666666666666666*delta*s*pow(tmax,3)*HeavisideTheta(-mrho + sqrt(s)))/pow(pow(mrho,3) - 1.*mrho*s,2) -
	       (0.666667*delta*pow(tmax,3)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,4) - 1.*pow(mrho,2)*s) +
	       (0.3333335*pow(delta,2)*s*pow(tmax,3)*HeavisideTheta(-mrho + sqrt(s)))/(pow(mrho,6) - 1.*pow(mrho,4)*s) -
	       0.0625*pow(eta1 - 1.*eta2,2)*(eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) +
	             pow(mpion,2)*(4.*pow(mpion,4) - 4.*pow(mrho,4) - 2.*pow(mpion,2)*s + 2.*pow(mrho,2)*s) +
	             pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
	          pow(eta1,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) + pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(pow(mrho,2) - 1.*s)*s + pow(mpion,4)*(-3.*pow(mrho,2) + s) +
	             pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
	             pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2))) +
	          pow(eta2,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) - 1.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(3.*pow(mrho,2) + s) +
	             pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
	             pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2))))*log(fabs(-1.*pow(ma1,2) + tmin)) +
	       (0.25*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-0.5*pow(ma1,6) - 0.5*pow(mpion,6) + 0.5*pow(mpion,4)*pow(mrho,2) +
	               pow(ma1,4)*(0.5*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) + pow(ma1,2)*pow(mpion,2)*(0.5*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s)) +
	            eta1*(pow(ma1,4)*(1.*pow(mpion,2) + 0.5*s) + pow(mpion,2)*(1.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-1.*pow(mrho,2) + 0.5*s) - 0.5*pow(mrho,2)*s) +
	               pow(ma1,2)*(-2.*pow(mpion,4) - 0.5*pow(mrho,2)*s + pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*log(fabs(-1.*pow(ma1,2) + tmin)))/(pow(ma1,2) - 1.*pow(mpion,2)) -
	       (0.25*(1.*eta1 - 1.*eta2)*(delta*(eta1*(1.*pow(ma1,6) - 1.*pow(mpion,6) + pow(mpion,4)*(-2.5*pow(mrho,2) + 0.5*s) + pow(mpion,2)*s*(-0.5*pow(mrho,2) + 1.*s) +
	                  pow(ma1,4)*(-3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s) + s*(0.5*pow(mrho,4) - 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
	                  pow(ma1,2)*(3.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s + 1.*pow(s,2))) +
	               eta2*(-1.*pow(ma1,6) + pow(ma1,4)*(3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
	                  pow(mpion,2)*(1.*pow(mpion,4) - 0.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
	                  pow(ma1,2)*(-3.*pow(mpion,4) + 1.*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.75*pow(s,2) + pow(mpion,2)*(1.*pow(mrho,2) + 1.*s)))) +
	            pow(mrho,2)*(eta2*(pow(ma1,4)*(-1. + 2.*C4*pow(mrho,2)) + pow(mpion,2)*(0.5*pow(mrho,2) + pow(mpion,2)*(-1. + 2.*C4*pow(mrho,2)) + 0.5*s) +
	                  pow(ma1,2)*(2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) + pow(mrho,2)*(-1.5 - 4.*C4*s) + s*(0.5 + 2.*C4*s))) +
	               eta1*(pow(ma1,4)*(1. - 2.*C4*pow(mrho,2)) + pow(mpion,4)*(1. - 2.*C4*pow(mrho,2)) + (-0.5*pow(mrho,2) + 0.5*s)*s +
	                  pow(ma1,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(-2. + 4.*C4*pow(mrho,2)) - 2.*C4*pow(s,2)) +
	                  pow(mpion,2)*(-4.*C4*pow(mrho,4) - 1.*s + pow(mrho,2)*(2. + 4.*C4*s)))))*log(fabs(-1.*pow(ma1,2) + tmin)))/pow(mrho,2) -
	       (HeavisideTheta(-mrho + sqrt(s))*(-0.5*(eta1*eta2*(0.25*pow(mrho,6) + 1.5*pow(mrho,4)*s - 0.125*delta*pow(mrho,4)*s - 0.75*pow(mrho,2)*pow(s,2) -
	                  0.75*delta*pow(mrho,2)*pow(s,2) + 0.375*delta*pow(s,3) + pow(ma1,4)*(-2.*pow(mrho,2) + 1.*delta*s) + pow(mpion,4)*(-6.*pow(mrho,2) + 3.*delta*s) +
	                  pow(mpion,2)*(-3.*pow(mrho,4) + (3. + 1.5*delta)*pow(mrho,2)*s - 1.5*delta*pow(s,2)) +
	                  pow(ma1,2)*(0.5*pow(mrho,4) + (-2.5 - 0.25*delta)*pow(mrho,2)*s + 1.25*delta*pow(s,2) + pow(mpion,2)*(6.*pow(mrho,2) - 3.*delta*s))) +
	               pow(eta1,2)*(0.5*pow(mrho,6) - 2.*pow(mrho,4)*s - 0.25*delta*pow(mrho,4)*s + 0.5*pow(mrho,2)*pow(s,2) + 1.*delta*pow(mrho,2)*pow(s,2) -
	                  0.25*delta*pow(s,3) + pow(mpion,4)*(3.*pow(mrho,2) - 1.5*delta*s) + pow(ma1,4)*(1.*pow(mrho,2) - 0.5*delta*s) +
	                  pow(mpion,2)*(4.*pow(mrho,4) + (-2. - 2.*delta)*pow(mrho,2)*s + 1.*delta*pow(s,2)) +
	                  pow(ma1,2)*(-1.5*pow(mrho,4) + (1.5 + 0.75*delta)*pow(mrho,2)*s - 0.75*delta*pow(s,2) + pow(mpion,2)*(-3.*pow(mrho,2) + 1.5*delta*s))) +
	               pow(eta2,2)*(-0.75*pow(mrho,6) + 0.5*pow(mrho,4)*s + 0.375*delta*pow(mrho,4)*s + 0.25*pow(mrho,2)*pow(s,2) - 0.25*delta*pow(mrho,2)*pow(s,2) -
	                  0.125*delta*pow(s,3) + pow(mpion,4)*(3.*pow(mrho,2) - 1.5*delta*s) + pow(ma1,4)*(1.*pow(mrho,2) - 0.5*delta*s) +
	                  pow(mpion,2)*(-1.*pow(mrho,4) + (-1. + 0.5*delta)*pow(mrho,2)*s + 0.5*delta*pow(s,2)) +
	                  pow(ma1,2)*(1.*pow(mrho,4) + (1. - 0.5*delta)*pow(mrho,2)*s - 0.5*delta*pow(s,2) + pow(mpion,2)*(-3.*pow(mrho,2) + 1.5*delta*s))))*tmin -
	            0.25*(eta1*eta2*(0.5*pow(mrho,4) - 2.5*pow(mrho,2)*s - 0.25*delta*pow(mrho,2)*s + 1.25*delta*pow(s,2) + pow(mpion,2)*(6.*pow(mrho,2) - 3.*delta*s) +
	                  pow(ma1,2)*(-2.*pow(mrho,2) + 1.*delta*s)) + pow(eta1,2)*
	                (-1.5*pow(mrho,4) + 1.5*pow(mrho,2)*s + 0.75*delta*pow(mrho,2)*s - 0.75*delta*pow(s,2) + pow(ma1,2)*(1.*pow(mrho,2) - 0.5*delta*s) +
	                  pow(mpion,2)*(-3.*pow(mrho,2) + 1.5*delta*s)) + pow(eta2,2)*
	                (1.*pow(mrho,4) + 1.*pow(mrho,2)*s - 0.5*delta*pow(mrho,2)*s - 0.5*delta*pow(s,2) + pow(ma1,2)*(1.*pow(mrho,2) - 0.5*delta*s) +
	                  pow(mpion,2)*(-3.*pow(mrho,2) + 1.5*delta*s)))*pow(tmin,2) -
	            0.16666666666666666*(pow(eta1,2)*(1.*pow(mrho,2) - 0.5*delta*s) + pow(eta2,2)*(1.*pow(mrho,2) - 0.5*delta*s) + eta1*eta2*(-2.*pow(mrho,2) + 1.*delta*s))*
	             pow(tmin,3) - 0.5*(eta1*eta2*(pow(mpion,6)*(2.*pow(mrho,2) - 1.*delta*s) + pow(ma1,6)*(-2.*pow(mrho,2) + 1.*delta*s) +
	                  pow(mpion,4)*(2.5*pow(mrho,4) + (-0.5 - 1.25*delta)*pow(mrho,2)*s + 0.25*delta*pow(s,2)) +
	                  s*(-0.25*pow(mrho,6) + 0.125*delta*pow(mrho,4)*s + 0.25*pow(mrho,2)*pow(s,2) - 0.125*delta*pow(s,3)) +
	                  pow(mpion,2)*(-0.25*pow(mrho,6) + (0.5 + 0.125*delta)*pow(mrho,4)*s + (-1.25 - 0.25*delta)*pow(mrho,2)*pow(s,2) + 0.625*delta*pow(s,3)) +
	                  pow(ma1,4)*(0.5*pow(mrho,4) + (-2.5 - 0.25*delta)*pow(mrho,2)*s + 1.25*delta*pow(s,2) + pow(mpion,2)*(6.*pow(mrho,2) - 3.*delta*s)) +
	                  pow(ma1,2)*(0.25*pow(mrho,6) + (1.5 - 0.125*delta)*pow(mrho,4)*s + (-0.75 - 0.75*delta)*pow(mrho,2)*pow(s,2) + 0.375*delta*pow(s,3) +
	                     pow(mpion,4)*(-6.*pow(mrho,2) + 3.*delta*s) + pow(mpion,2)*(-3.*pow(mrho,4) + (3. + 1.5*delta)*pow(mrho,2)*s - 1.5*delta*pow(s,2)))) +
	               pow(eta2,2)*(pow(ma1,6)*(1.*pow(mrho,2) - 0.5*delta*s) +
	                  pow(mpion,2)*(0.25*pow(mrho,6) + (-0.5 - 0.125*delta)*pow(mrho,4)*s + (0.25 + 0.25*delta)*pow(mrho,2)*pow(s,2) - 0.125*delta*pow(s,3) +
	                     pow(mpion,4)*(-1.*pow(mrho,2) + 0.5*delta*s)) + pow(ma1,4)*
	                   (1.*pow(mrho,4) + (1. - 0.5*delta)*pow(mrho,2)*s - 0.5*delta*pow(s,2) + pow(mpion,2)*(-3.*pow(mrho,2) + 1.5*delta*s)) +
	                  pow(ma1,2)*(-0.75*pow(mrho,6) + (0.5 + 0.375*delta)*pow(mrho,4)*s + (0.25 - 0.25*delta)*pow(mrho,2)*pow(s,2) - 0.125*delta*pow(s,3) +
	                     pow(mpion,4)*(3.*pow(mrho,2) - 1.5*delta*s) + pow(mpion,2)*(-1.*pow(mrho,4) + (-1. + 0.5*delta)*pow(mrho,2)*s + 0.5*delta*pow(s,2)))) +
	               pow(eta1,2)*(pow(ma1,6)*(1.*pow(mrho,2) - 0.5*delta*s) + pow(mpion,2)*pow(s,2)*(1.*pow(mrho,2) - 0.5*delta*s) +
	                  pow(mpion,6)*(-1.*pow(mrho,2) + 0.5*delta*s) + pow(mpion,4)*(-2.5*pow(mrho,4) + (0.5 + 1.25*delta)*pow(mrho,2)*s - 0.25*delta*pow(s,2)) +
	                  s*(0.25*pow(mrho,6) - 0.125*delta*pow(mrho,4)*s - 0.25*pow(mrho,2)*pow(s,2) + 0.125*delta*pow(s,3)) +
	                  pow(ma1,4)*(-1.5*pow(mrho,4) + (1.5 + 0.75*delta)*pow(mrho,2)*s - 0.75*delta*pow(s,2) + pow(mpion,2)*(-3.*pow(mrho,2) + 1.5*delta*s)) +
	                  pow(ma1,2)*(0.5*pow(mrho,6) + (-2. - 0.25*delta)*pow(mrho,4)*s + (0.5 + 1.*delta)*pow(mrho,2)*pow(s,2) - 0.25*delta*pow(s,3) +
	                     pow(mpion,4)*(3.*pow(mrho,2) - 1.5*delta*s) + pow(mpion,2)*(4.*pow(mrho,4) + (-2. - 2.*delta)*pow(mrho,2)*s + 1.*delta*pow(s,2)))))*
	             log(fabs(-1.*pow(ma1,2) + tmin))))/(pow(mrho,2)*(pow(mrho,2) - 1.*s)) - 0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-1.*pow(mpion,2) + tmin)) +
	       (0.5*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*pow(mpion,4)*(-1.*pow(mrho,2) + 1.*s) +
	            eta1*pow(mpion,2)*(-0.5*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + 0.5*pow(mrho,2)*s))*log(fabs(-1.*pow(mpion,2) + tmin)))/(pow(ma1,2) - 1.*pow(mpion,2)) -
	       (0.5*(pow(mrho,2) - 1.*s)*((0.5 - 0.25*delta)*pow(mrho,2) + C4*(-2. + 1.*delta)*pow(mrho,4) + (-0.25 + 0.125*delta)*delta*s)*pow(tmin,2) +
	          tmin*(-0.5*pow(mrho,6) + 0.25*delta*pow(mrho,6) + 2.*C4*pow(mrho,8) - 1.*C4*delta*pow(mrho,8) + 1.*pow(mrho,4)*s + 0.25*delta*pow(mrho,4)*s -
	             0.375*pow(delta,2)*pow(mrho,4)*s - 6.*C4*pow(mrho,6)*s + 3.*C4*delta*pow(mrho,6)*s - 0.5*pow(mrho,2)*pow(s,2) - 0.75*delta*pow(mrho,2)*pow(s,2) +
	             0.5*pow(delta,2)*pow(mrho,2)*pow(s,2) + 4.*C4*pow(mrho,4)*pow(s,2) - 2.*C4*delta*pow(mrho,4)*pow(s,2) + 0.25*delta*pow(s,3) -
	             0.125*pow(delta,2)*pow(s,3) + pow(mpion,2)*(C4*(2. - 1.*delta)*pow(mrho,6) + (0.5 - 0.125*pow(delta,2))*pow(mrho,2)*s + (-1.25 + 0.625*delta)*delta*pow(s,2) +
	                pow(mrho,4)*(-0.5 + 1.25*delta - 0.5*pow(delta,2) - 2.*C4*s + 1.*C4*delta*s)) +
	             (pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) - 0.5000000000000001*pow(2. - 1.*delta,2)*pow(mrho,2)*s + (1. - 0.5*delta)*delta*pow(s,2)) +
	                s*((-0.5 + 0.25*delta)*pow(mrho,4) + (0.5 - 0.125*pow(delta,2))*pow(mrho,2)*s + (-0.25 + 0.125*delta)*delta*pow(s,2)))*HeavisideTheta(-mrho + sqrt(s))) +
	          pow(mrho,2)*(s*((0.5 - 0.75*delta + 0.25*pow(delta,2))*pow(mrho,4) - 0.12500000000000003*pow(2. - 1.*delta,2)*pow(mrho,2)*s +
	                (0.25 - 0.125*delta)*delta*pow(s,2)) + pow(mpion,2)*(C4*(4. - 2.*delta)*pow(mrho,6) + delta*(-3. + 1.5*delta)*pow(s,2) +
	                pow(mrho,2)*s*(2. - 1.5*pow(delta,2) + 4.*C4*s + delta*(2. - 2.*C4*s)) + pow(mrho,4)*(-2. - 8.*C4*s + delta*(1. + 4.*C4*s))) +
	             s*((0.5 - 0.25*delta)*pow(mrho,4) + 0.12500000000000003*pow(2. - 1.*delta,2)*pow(mrho,2)*s + (-0.25 + 0.125*delta)*delta*pow(s,2) +
	                pow(mpion,2)*((-4. + 2.*delta)*pow(mrho,2) + (2. - 1.*delta)*delta*s))*HeavisideTheta(-mrho + sqrt(s)))*log(fabs(-1.*pow(mpion,2) + tmin)))/
	        (pow(mrho,4)*(pow(mrho,2) - 1.*s)) + 0.0625*pow(eta1 - 1.*eta2,2)*
	        (eta1*eta2*(-4.*pow(ma1,6) + pow(ma1,4)*(12.*pow(mpion,2) - 6.*s) + pow(mpion,2)*(4.*pow(mpion,4) - 4.*pow(mrho,4) - 2.*pow(mpion,2)*s + 2.*pow(mrho,2)*s) +
	             pow(ma1,2)*(-12.*pow(mpion,4) + 2.*pow(mrho,4) + 8.*pow(mpion,2)*s + 4.*pow(mrho,2)*s - 4.*pow(s,2))) +
	          pow(eta1,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) + pow(mpion,2)*pow(mrho,2)*s + pow(mrho,2)*(pow(mrho,2) - 1.*s)*s + pow(mpion,4)*(-3.*pow(mrho,2) + s) +
	             pow(ma1,4)*(-6.*pow(mpion,2) - 3.*pow(mrho,2) + 3.*s) +
	             pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(6.*pow(mrho,2) - 4.*s) - 4.*pow(mrho,2)*s + 2.*pow(s,2))) +
	          pow(eta2,2)*(2.*pow(ma1,6) - 2.*pow(mpion,6) - 1.*pow(mpion,2)*pow(mrho,2)*s + pow(mpion,4)*(3.*pow(mrho,2) + s) +
	             pow(ma1,4)*(-6.*pow(mpion,2) + 3.*pow(mrho,2) + 3.*s) +
	             pow(ma1,2)*(6.*pow(mpion,4) + pow(mrho,4) + pow(mpion,2)*(-6.*pow(mrho,2) - 4.*s) - 2.*pow(mrho,2)*s + 2.*pow(s,2))))*log(fabs(-1.*pow(ma1,2) + tmax)) -
	       (0.25*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*(-0.5*pow(ma1,6) - 0.5*pow(mpion,6) + 0.5*pow(mpion,4)*pow(mrho,2) +
	               pow(ma1,4)*(0.5*pow(mpion,2) + 0.5*pow(mrho,2) - 1.*s) + pow(ma1,2)*pow(mpion,2)*(0.5*pow(mpion,2) + 1.*pow(mrho,2) - 1.*s)) +
	            eta1*(pow(ma1,4)*(1.*pow(mpion,2) + 0.5*s) + pow(mpion,2)*(1.*pow(mpion,4) + 1.*pow(mrho,4) + pow(mpion,2)*(-1.*pow(mrho,2) + 0.5*s) - 0.5*pow(mrho,2)*s) +
	               pow(ma1,2)*(-2.*pow(mpion,4) - 0.5*pow(mrho,2)*s + pow(mpion,2)*(-1.*pow(mrho,2) + 1.*s))))*log(fabs(-1.*pow(ma1,2) + tmax)))/(pow(ma1,2) - 1.*pow(mpion,2)) +
	       (0.25*(1.*eta1 - 1.*eta2)*(delta*(eta1*(1.*pow(ma1,6) - 1.*pow(mpion,6) + pow(mpion,4)*(-2.5*pow(mrho,2) + 0.5*s) + pow(mpion,2)*s*(-0.5*pow(mrho,2) + 1.*s) +
	                  pow(ma1,4)*(-3.*pow(mpion,2) - 1.5*pow(mrho,2) + 1.5*s) + s*(0.5*pow(mrho,4) - 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
	                  pow(ma1,2)*(3.*pow(mpion,4) + 0.5*pow(mrho,4) + pow(mpion,2)*(4.*pow(mrho,2) - 2.*s) - 2.*pow(mrho,2)*s + 1.*pow(s,2))) +
	               eta2*(-1.*pow(ma1,6) + pow(ma1,4)*(3.*pow(mpion,2) - 1.*pow(mrho,2) - 1.*s) +
	                  pow(mpion,2)*(1.*pow(mpion,4) - 0.5*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.25*pow(s,2)) +
	                  pow(ma1,2)*(-3.*pow(mpion,4) + 1.*pow(mrho,4) + 0.25*pow(mrho,2)*s - 0.75*pow(s,2) + pow(mpion,2)*(1.*pow(mrho,2) + 1.*s)))) +
	            pow(mrho,2)*(eta2*(pow(ma1,4)*(-1. + 2.*C4*pow(mrho,2)) + pow(mpion,2)*(0.5*pow(mrho,2) + pow(mpion,2)*(-1. + 2.*C4*pow(mrho,2)) + 0.5*s) +
	                  pow(ma1,2)*(2.*C4*pow(mrho,4) + pow(mpion,2)*(2. - 4.*C4*pow(mrho,2)) + pow(mrho,2)*(-1.5 - 4.*C4*s) + s*(0.5 + 2.*C4*s))) +
	               eta1*(pow(ma1,4)*(1. - 2.*C4*pow(mrho,2)) + pow(mpion,4)*(1. - 2.*C4*pow(mrho,2)) + (-0.5*pow(mrho,2) + 0.5*s)*s +
	                  pow(ma1,2)*(-1.*pow(mrho,2) + 2.*C4*pow(mrho,4) + pow(mpion,2)*(-2. + 4.*C4*pow(mrho,2)) - 2.*C4*pow(s,2)) +
	                  pow(mpion,2)*(-4.*C4*pow(mrho,4) - 1.*s + pow(mrho,2)*(2. + 4.*C4*s)))))*log(fabs(-1.*pow(ma1,2) + tmax)))/pow(mrho,2) +
	       (HeavisideTheta(-mrho + sqrt(s))*(-0.5*(eta1*eta2*(0.25*pow(mrho,6) + 1.5*pow(mrho,4)*s - 0.125*delta*pow(mrho,4)*s - 0.75*pow(mrho,2)*pow(s,2) -
	                  0.75*delta*pow(mrho,2)*pow(s,2) + 0.375*delta*pow(s,3) + pow(ma1,4)*(-2.*pow(mrho,2) + 1.*delta*s) + pow(mpion,4)*(-6.*pow(mrho,2) + 3.*delta*s) +
	                  pow(mpion,2)*(-3.*pow(mrho,4) + (3. + 1.5*delta)*pow(mrho,2)*s - 1.5*delta*pow(s,2)) +
	                  pow(ma1,2)*(0.5*pow(mrho,4) + (-2.5 - 0.25*delta)*pow(mrho,2)*s + 1.25*delta*pow(s,2) + pow(mpion,2)*(6.*pow(mrho,2) - 3.*delta*s))) +
	               pow(eta1,2)*(0.5*pow(mrho,6) - 2.*pow(mrho,4)*s - 0.25*delta*pow(mrho,4)*s + 0.5*pow(mrho,2)*pow(s,2) + 1.*delta*pow(mrho,2)*pow(s,2) -
	                  0.25*delta*pow(s,3) + pow(mpion,4)*(3.*pow(mrho,2) - 1.5*delta*s) + pow(ma1,4)*(1.*pow(mrho,2) - 0.5*delta*s) +
	                  pow(mpion,2)*(4.*pow(mrho,4) + (-2. - 2.*delta)*pow(mrho,2)*s + 1.*delta*pow(s,2)) +
	                  pow(ma1,2)*(-1.5*pow(mrho,4) + (1.5 + 0.75*delta)*pow(mrho,2)*s - 0.75*delta*pow(s,2) + pow(mpion,2)*(-3.*pow(mrho,2) + 1.5*delta*s))) +
	               pow(eta2,2)*(-0.75*pow(mrho,6) + 0.5*pow(mrho,4)*s + 0.375*delta*pow(mrho,4)*s + 0.25*pow(mrho,2)*pow(s,2) - 0.25*delta*pow(mrho,2)*pow(s,2) -
	                  0.125*delta*pow(s,3) + pow(mpion,4)*(3.*pow(mrho,2) - 1.5*delta*s) + pow(ma1,4)*(1.*pow(mrho,2) - 0.5*delta*s) +
	                  pow(mpion,2)*(-1.*pow(mrho,4) + (-1. + 0.5*delta)*pow(mrho,2)*s + 0.5*delta*pow(s,2)) +
	                  pow(ma1,2)*(1.*pow(mrho,4) + (1. - 0.5*delta)*pow(mrho,2)*s - 0.5*delta*pow(s,2) + pow(mpion,2)*(-3.*pow(mrho,2) + 1.5*delta*s))))*tmax -
	            0.25*(eta1*eta2*(0.5*pow(mrho,4) - 2.5*pow(mrho,2)*s - 0.25*delta*pow(mrho,2)*s + 1.25*delta*pow(s,2) + pow(mpion,2)*(6.*pow(mrho,2) - 3.*delta*s) +
	                  pow(ma1,2)*(-2.*pow(mrho,2) + 1.*delta*s)) + pow(eta1,2)*
	                (-1.5*pow(mrho,4) + 1.5*pow(mrho,2)*s + 0.75*delta*pow(mrho,2)*s - 0.75*delta*pow(s,2) + pow(ma1,2)*(1.*pow(mrho,2) - 0.5*delta*s) +
	                  pow(mpion,2)*(-3.*pow(mrho,2) + 1.5*delta*s)) + pow(eta2,2)*
	                (1.*pow(mrho,4) + 1.*pow(mrho,2)*s - 0.5*delta*pow(mrho,2)*s - 0.5*delta*pow(s,2) + pow(ma1,2)*(1.*pow(mrho,2) - 0.5*delta*s) +
	                  pow(mpion,2)*(-3.*pow(mrho,2) + 1.5*delta*s)))*pow(tmax,2) -
	            0.16666666666666666*(pow(eta1,2)*(1.*pow(mrho,2) - 0.5*delta*s) + pow(eta2,2)*(1.*pow(mrho,2) - 0.5*delta*s) + eta1*eta2*(-2.*pow(mrho,2) + 1.*delta*s))*
	             pow(tmax,3) - 0.5*(eta1*eta2*(pow(mpion,6)*(2.*pow(mrho,2) - 1.*delta*s) + pow(ma1,6)*(-2.*pow(mrho,2) + 1.*delta*s) +
	                  pow(mpion,4)*(2.5*pow(mrho,4) + (-0.5 - 1.25*delta)*pow(mrho,2)*s + 0.25*delta*pow(s,2)) +
	                  s*(-0.25*pow(mrho,6) + 0.125*delta*pow(mrho,4)*s + 0.25*pow(mrho,2)*pow(s,2) - 0.125*delta*pow(s,3)) +
	                  pow(mpion,2)*(-0.25*pow(mrho,6) + (0.5 + 0.125*delta)*pow(mrho,4)*s + (-1.25 - 0.25*delta)*pow(mrho,2)*pow(s,2) + 0.625*delta*pow(s,3)) +
	                  pow(ma1,4)*(0.5*pow(mrho,4) + (-2.5 - 0.25*delta)*pow(mrho,2)*s + 1.25*delta*pow(s,2) + pow(mpion,2)*(6.*pow(mrho,2) - 3.*delta*s)) +
	                  pow(ma1,2)*(0.25*pow(mrho,6) + (1.5 - 0.125*delta)*pow(mrho,4)*s + (-0.75 - 0.75*delta)*pow(mrho,2)*pow(s,2) + 0.375*delta*pow(s,3) +
	                     pow(mpion,4)*(-6.*pow(mrho,2) + 3.*delta*s) + pow(mpion,2)*(-3.*pow(mrho,4) + (3. + 1.5*delta)*pow(mrho,2)*s - 1.5*delta*pow(s,2)))) +
	               pow(eta2,2)*(pow(ma1,6)*(1.*pow(mrho,2) - 0.5*delta*s) +
	                  pow(mpion,2)*(0.25*pow(mrho,6) + (-0.5 - 0.125*delta)*pow(mrho,4)*s + (0.25 + 0.25*delta)*pow(mrho,2)*pow(s,2) - 0.125*delta*pow(s,3) +
	                     pow(mpion,4)*(-1.*pow(mrho,2) + 0.5*delta*s)) + pow(ma1,4)*
	                   (1.*pow(mrho,4) + (1. - 0.5*delta)*pow(mrho,2)*s - 0.5*delta*pow(s,2) + pow(mpion,2)*(-3.*pow(mrho,2) + 1.5*delta*s)) +
	                  pow(ma1,2)*(-0.75*pow(mrho,6) + (0.5 + 0.375*delta)*pow(mrho,4)*s + (0.25 - 0.25*delta)*pow(mrho,2)*pow(s,2) - 0.125*delta*pow(s,3) +
	                     pow(mpion,4)*(3.*pow(mrho,2) - 1.5*delta*s) + pow(mpion,2)*(-1.*pow(mrho,4) + (-1. + 0.5*delta)*pow(mrho,2)*s + 0.5*delta*pow(s,2)))) +
	               pow(eta1,2)*(pow(ma1,6)*(1.*pow(mrho,2) - 0.5*delta*s) + pow(mpion,2)*pow(s,2)*(1.*pow(mrho,2) - 0.5*delta*s) +
	                  pow(mpion,6)*(-1.*pow(mrho,2) + 0.5*delta*s) + pow(mpion,4)*(-2.5*pow(mrho,4) + (0.5 + 1.25*delta)*pow(mrho,2)*s - 0.25*delta*pow(s,2)) +
	                  s*(0.25*pow(mrho,6) - 0.125*delta*pow(mrho,4)*s - 0.25*pow(mrho,2)*pow(s,2) + 0.125*delta*pow(s,3)) +
	                  pow(ma1,4)*(-1.5*pow(mrho,4) + (1.5 + 0.75*delta)*pow(mrho,2)*s - 0.75*delta*pow(s,2) + pow(mpion,2)*(-3.*pow(mrho,2) + 1.5*delta*s)) +
	                  pow(ma1,2)*(0.5*pow(mrho,6) + (-2. - 0.25*delta)*pow(mrho,4)*s + (0.5 + 1.*delta)*pow(mrho,2)*pow(s,2) - 0.25*delta*pow(s,3) +
	                     pow(mpion,4)*(3.*pow(mrho,2) - 1.5*delta*s) + pow(mpion,2)*(4.*pow(mrho,4) + (-2. - 2.*delta)*pow(mrho,2)*s + 1.*delta*pow(s,2)))))*
	             log(fabs(-1.*pow(ma1,2) + tmax))))/(pow(mrho,2)*(pow(mrho,2) - 1.*s)) + 0.5*pow(-2. + delta,2)*pow(mpion,2)*log(fabs(-1.*pow(mpion,2) + tmax)) -
	       (0.5*(-2. + delta)*(eta1 - 1.*eta2)*(eta2*pow(mpion,4)*(-1.*pow(mrho,2) + 1.*s) +
	            eta1*pow(mpion,2)*(-0.5*pow(mrho,4) + pow(mpion,2)*(1.*pow(mrho,2) - 1.*s) + 0.5*pow(mrho,2)*s))*log(fabs(-1.*pow(mpion,2) + tmax)))/(pow(ma1,2) - 1.*pow(mpion,2)) +
	       (0.5*(pow(mrho,2) - 1.*s)*((0.5 - 0.25*delta)*pow(mrho,2) + C4*(-2. + 1.*delta)*pow(mrho,4) + (-0.25 + 0.125*delta)*delta*s)*pow(tmax,2) +
	          tmax*(-0.5*pow(mrho,6) + 0.25*delta*pow(mrho,6) + 2.*C4*pow(mrho,8) - 1.*C4*delta*pow(mrho,8) + 1.*pow(mrho,4)*s + 0.25*delta*pow(mrho,4)*s -
	             0.375*pow(delta,2)*pow(mrho,4)*s - 6.*C4*pow(mrho,6)*s + 3.*C4*delta*pow(mrho,6)*s - 0.5*pow(mrho,2)*pow(s,2) - 0.75*delta*pow(mrho,2)*pow(s,2) +
	             0.5*pow(delta,2)*pow(mrho,2)*pow(s,2) + 4.*C4*pow(mrho,4)*pow(s,2) - 2.*C4*delta*pow(mrho,4)*pow(s,2) + 0.25*delta*pow(s,3) -
	             0.125*pow(delta,2)*pow(s,3) + pow(mpion,2)*(C4*(2. - 1.*delta)*pow(mrho,6) + (0.5 - 0.125*pow(delta,2))*pow(mrho,2)*s + (-1.25 + 0.625*delta)*delta*pow(s,2) +
	                pow(mrho,4)*(-0.5 + 1.25*delta - 0.5*pow(delta,2) - 2.*C4*s + 1.*C4*delta*s)) +
	             (pow(mpion,2)*((-2. + 1.*delta)*pow(mrho,4) - 0.5000000000000001*pow(2. - 1.*delta,2)*pow(mrho,2)*s + (1. - 0.5*delta)*delta*pow(s,2)) +
	                s*((-0.5 + 0.25*delta)*pow(mrho,4) + (0.5 - 0.125*pow(delta,2))*pow(mrho,2)*s + (-0.25 + 0.125*delta)*delta*pow(s,2)))*HeavisideTheta(-mrho + sqrt(s))) +
	          pow(mrho,2)*(s*((0.5 - 0.75*delta + 0.25*pow(delta,2))*pow(mrho,4) - 0.12500000000000003*pow(2. - 1.*delta,2)*pow(mrho,2)*s +
	                (0.25 - 0.125*delta)*delta*pow(s,2)) + pow(mpion,2)*(C4*(4. - 2.*delta)*pow(mrho,6) + delta*(-3. + 1.5*delta)*pow(s,2) +
	                pow(mrho,2)*s*(2. - 1.5*pow(delta,2) + 4.*C4*s + delta*(2. - 2.*C4*s)) + pow(mrho,4)*(-2. - 8.*C4*s + delta*(1. + 4.*C4*s))) +
	             s*((0.5 - 0.25*delta)*pow(mrho,4) + 0.12500000000000003*pow(2. - 1.*delta,2)*pow(mrho,2)*s + (-0.25 + 0.125*delta)*delta*pow(s,2) +
	                pow(mpion,2)*((-4. + 2.*delta)*pow(mrho,2) + (2. - 1.*delta)*delta*s))*HeavisideTheta(-mrho + sqrt(s)))*log(fabs(-1.*pow(mpion,2) + tmax)))/
	        (pow(mrho,4)*(pow(mrho,2) - 1.*s))))/(16.*Pi*s*(-4*pow(mpion,2) + s));

  return sigma / spin_deg_factor;
}
