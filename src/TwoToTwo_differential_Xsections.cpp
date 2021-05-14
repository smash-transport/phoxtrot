/*
 * Collection of differential cross sections as a function of mandelstam t, the
 * center-of-mass energy (s). The cross section does in principle also depend
 * on the mass of the rho meson, it's value is defined in "common_functions.h".
 * The naming conventions (C11, C12, ...) follow from the PhD thesis of
 * S. Turbide (See README for link).
 */

#include "include/TwoToTwo_differential_Xsections.h"
#include "include/common_functions.h"

// C11
double TwoToTwo_Diff_Xsections::diff_xsection_C11(double t, double s) {
  double diff_sigma;
  double spin_deg_factor = 3.0;

  diff_sigma =
      (pow(Const, 2) * pow(ghat, 4) *
       ((-8 * pow(-2 + delta, 2) * pow(mpion, 2)) /
            (pow(mrho, 2) * pow(pow(mpion, 2) - s, 2)) -
        (8 * pow(-2 + delta, 2) * pow(mpion, 2) *
         (pow(mpion, 4) + pow(pow(mrho, 2) - t, 2) -
          2 * pow(mpion, 2) * (pow(mrho, 2) + t))) /
            (pow(mrho, 2) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             pow(pow(mpion, 2) - t, 2)) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (-(eta2 * (pow(mpion, 2) + s)) + eta1 * (-pow(mrho, 2) + s + t)) *
         (-pow(mpion, 4) + pow(mpion, 2) * (pow(mrho, 2) - 2 * t) +
          t * (-pow(mrho, 2) + 2 * s + t))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(mpion, 2) - t)) -
        (8 * (-2 + delta) *
         (pow(mpion, 4) * (2 - 3 * delta + 8 * C4 * pow(mrho, 2)) +
          pow(mrho, 4) * (-2 + delta + 8 * C4 * t) +
          t * ((2 + 3 * delta) * s + 2 * delta * t) +
          pow(mpion, 2) *
              (-8 * C4 * pow(mrho, 4) + (-2 + delta) * s - (2 + 3 * delta) * t +
               4 * pow(mrho, 2) * (1 + 4 * C4 * t)) -
          pow(mrho, 2) * (t * (-2 + 3 * delta + 8 * C4 * t) +
                          s * (-2 + delta + 16 * C4 * t)))) /
            (pow(mrho, 2) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(mpion, 2) - t)) +
        (4 * (-2 + delta) * (eta1 - eta2) * (pow(ma1, 2) - s) *
         (eta2 * (pow(mpion, 2) + s) *
              (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
               s * (pow(mrho, 2) - s - 2 * t)) +
          eta1 * (-4 * pow(mpion, 6) +
                  s * (-pow(mrho, 2) + s) * (-pow(mrho, 2) + s + t) +
                  pow(mpion, 4) * (3 * pow(mrho, 2) + s + t) -
                  pow(mpion, 2) * (pow(mrho, 4) + 2 * s * (s - t) +
                                   pow(mrho, 2) * (-s + t))))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (pow(eta2, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 4) * (pow(pow(mrho, 2) + 2 * s, 2) - 2 * s * t) +
               pow(s, 2) * (pow(pow(mrho, 2) + s, 2) +
                            2 * (-pow(mrho, 2) + s) * t + 2 * pow(t, 2)) -
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * (2 * s - t) +
                    2 * s * (s + t))) +
          2 * eta1 * eta2 *
              (-pow(mpion, 8) +
               pow(mpion, 4) * (pow(mrho, 4) + 2 * pow(mrho, 2) * s +
                                2 * s * (-2 * s + t)) -
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * (s + t) - 2 * s * (s + t)) +
               pow(s, 2) * (pow(mrho, 4) - pow(s, 2) + 2 * pow(mrho, 2) * t -
                            2 * t * (s + t))) +
          pow(eta1, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 4) * (3 * pow(mrho, 4) + 2 * s * (2 * s - t) +
                                2 * pow(mrho, 2) * (-3 * s + t)) -
               2 * pow(mpion, 2) * (pow(mrho, 2) - s) *
                   (-2 * s * (s + t) + pow(mrho, 2) * (2 * s + t)) +
               s * (-pow(mrho, 2) + s) *
                   (pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                    pow(mrho, 2) * (s + 2 * t))))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(mpion, 8) -
               pow(mpion, 4) *
                   (pow(mrho, 4) + 2 * (pow(mrho, 2) + s) * t - 4 * pow(t, 2)) +
               pow(t, 2) * (-pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                            2 * pow(s, 2) + 2 * s * t + pow(t, 2)) +
               2 * pow(mpion, 2) * t *
                   (pow(mrho, 4) + pow(mrho, 2) * (s + t) - 2 * t * (s + t))) +
          pow(eta2, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 4) *
                   (pow(mrho, 4) + 4 * pow(mrho, 2) * t - 2 * (s - 2 * t) * t) +
               pow(t, 2) * (pow(mrho, 4) + 2 * pow(s, 2) + 2 * s * t +
                            pow(t, 2) + 2 * pow(mrho, 2) * (-s + t)) -
               2 * pow(mpion, 2) * t *
                   (pow(mrho, 4) - pow(mrho, 2) * (s - 2 * t) +
                    2 * t * (s + t))) +
          pow(eta1, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 4) *
                   (3 * pow(mrho, 4) + 2 * pow(mrho, 2) * (s - 3 * t) -
                    2 * (s - 2 * t) * t) +
               t * (-pow(mrho, 2) + t) *
                   (2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                    pow(mrho, 2) * (2 * s + t)) -
               2 * pow(mpion, 2) * (-pow(mrho, 2) + t) *
                   (2 * t * (s + t) - pow(mrho, 2) * (s + 2 * t))))) /
            ((pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             pow(pow(ma1, 2) - t, 2)) +
        (8 * (-2 + delta) *
         ((-2 + delta) * pow(mrho, 6) +
          pow(mpion, 6) * (-2 + 3 * delta - 8 * C4 * pow(mrho, 2)) +
          s * t * ((-2 + 3 * delta) * s + 4 * delta * t) +
          pow(mpion, 4) *
              (8 * C4 * pow(mrho, 4) + 4 * delta * s + 2 * t - 3 * delta * t -
               pow(mrho, 2) * (2 + delta + 16 * C4 * s - 8 * C4 * t)) +
          pow(mrho, 4) *
              (-((-2 + delta) * t) + s * (4 - 2 * delta + 8 * C4 * t)) +
          pow(mrho, 2) * s *
              (s * (-2 + delta - 8 * C4 * t) - 2 * t * (delta + 8 * C4 * t)) +
          pow(mpion, 2) *
              (s * ((2 - 3 * delta) * s - 8 * delta * t) -
               pow(mrho, 4) * (-6 + 3 * delta + 8 * C4 * (s + t)) +
               pow(mrho, 2) * (8 * C4 * pow(s, 2) + 4 * (-1 + delta) * t +
                               s * (-8 + 6 * delta + 32 * C4 * t))))) /
            (pow(mrho, 2) * (pow(mpion, 2) - s) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(mpion, 2) - t)) +
        (2 * pow(eta1 - eta2, 2) * (pow(ma1, 2) - s) *
         (pow(eta1, 2) * (pow(mpion, 8) +
                          pow(mpion, 4) * (2 * pow(mrho, 4) + 2 * s * t -
                                           3 * pow(mrho, 2) * (s + t)) +
                          s * t *
                              (2 * pow(mrho, 4) + pow(s, 2) + 3 * s * t +
                               pow(t, 2) - 3 * pow(mrho, 2) * (s + t)) -
                          2 * pow(mpion, 2) * (pow(mrho, 2) - s - t) *
                              (-2 * s * t + pow(mrho, 2) * (s + t))) +
          pow(eta2, 2) * (pow(mpion, 8) -
                          4 * pow(mpion, 2) * s * t * (pow(mrho, 2) + s + t) +
                          pow(mpion, 4) * (2 * s * t + pow(mrho, 2) * (s + t)) +
                          s * t *
                              (pow(s, 2) + 3 * s * t + pow(t, 2) +
                               pow(mrho, 2) * (s + t))) +
          2 * eta1 * eta2 *
              (-pow(mpion, 8) + 2 * pow(mpion, 6) * pow(mrho, 2) -
               2 * pow(mpion, 4) * s * t -
               s * t *
                   (pow(s, 2) + 3 * s * t + pow(t, 2) -
                    2 * pow(mrho, 2) * (s + t)) -
               pow(mpion, 2) *
                   (-4 * s * t * (s + t) +
                    pow(mrho, 2) * (pow(s, 2) + 4 * s * t + pow(t, 2)))))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
             (pow(ma1, 2) - t)) +
        (8 *
         (pow(delta, 2) * (8 * pow(mpion, 4) + 3 * pow(mrho, 4) -
                           6 * pow(mrho, 2) * (s + t) + 2 * pow(s + t, 2) +
                           4 * pow(mpion, 2) *
                               (3 * pow(mrho, 2) - 2 * (s + t))) -
          4 * delta * pow(mrho, 2) *
              (16 * C4 * pow(mpion, 4) + pow(mrho, 2) * (3 - 6 * C4 * (s + t)) +
               (s + t) * (-3 + 4 * C4 * (s + t)) +
               2 * pow(mpion, 2) *
                   (3 + C4 * (6 * pow(mrho, 2) - 8 * (s + t)))) +
          4 * pow(mrho, 4) *
              (3 + 4 * C4 * (2 * pow(mpion, 2) - s - t) *
                       (3 + C4 * (4 * pow(mpion, 2) - 2 * (s + t)))))) /
            (pow(mrho, 4) * (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
                             2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) * (-pow(ma1, 2) + s) *
         (eta2 * (-2 * pow(mpion, 4) * (delta - 4 * C4 * pow(mrho, 2)) *
                      (pow(mrho, 2) + 4 * s) +
                  pow(mpion, 2) *
                      (-2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * s) +
                       8 * delta * s * (s + t) -
                       pow(mrho, 2) * ((-10 + delta) * s - (-2 + delta) * t +
                                       32 * C4 * s * (s + t))) +
                  s * (2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * s) -
                       2 * delta * pow(s + t, 2) +
                       pow(mrho, 2) * ((-6 + delta) * s + (-2 + delta) * t +
                                       8 * C4 * pow(s + t, 2)))) +
          eta1 *
              (4 * pow(mpion, 4) *
                   (6 * C4 * pow(mrho, 4) + 2 * delta * s +
                    pow(mrho, 2) * (1 - 2 * delta - 8 * C4 * s)) +
               2 * delta * s * pow(s + t, 2) -
               pow(mrho, 2) *
                   ((-6 + 5 * delta) * pow(s, 2) +
                    2 * (-2 + 3 * delta) * s * t + (-2 + delta) * pow(t, 2) +
                    8 * C4 * s * pow(s + t, 2)) +
               pow(mrho, 4) *
                   ((-2 + delta) * (3 * s + t) + 8 * C4 * s * (s + 2 * t)) -
               2 * pow(mpion, 2) *
                   (4 * delta * s * (s + t) -
                    pow(mrho, 2) * (-6 * s + 7 * delta * s - 2 * t +
                                    3 * delta * t + 16 * C4 * s * (s + t)) +
                    2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * (2 * s + t)))))) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
              2 * pow(mpion, 2) * (pow(mrho, 2) + s))) +
        (4 * (eta1 - eta2) *
         (((-2 + delta) *
           (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
            s * (pow(mrho, 2) - s - 2 * t)) *
           (eta1 * (pow(mrho, 2) - s - t) + eta2 * (pow(mpion, 2) + t))) /
              ((pow(mpion, 2) - s) * (pow(ma1, 2) - t)) +
          ((-2 + delta) *
           (eta2 * (pow(mpion, 2) + t) *
                (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * t) +
                 (pow(mrho, 2) - 2 * s - t) * t) +
            eta1 * (-4 * pow(mpion, 6) +
                    (pow(mrho, 2) - t) * (pow(mrho, 2) - s - t) * t +
                    pow(mpion, 4) * (3 * pow(mrho, 2) + s + t) -
                    pow(mpion, 2) * (pow(mrho, 4) + pow(mrho, 2) * (s - t) +
                                     2 * t * (-s + t))))) /
              ((-pow(ma1, 2) + t) * (-pow(mpion, 2) + t)) +
          (eta2 *
               (-2 * pow(mpion, 4) * (delta - 4 * C4 * pow(mrho, 2)) *
                    (pow(mrho, 2) + 4 * t) +
                pow(mpion, 2) * (8 * delta * t * (s + t) -
                                 2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * t) -
                                 pow(mrho, 2) *
                                     (-((-2 + delta) * s) + (-10 + delta) * t +
                                      32 * C4 * t * (s + t))) +
                t * (-2 * delta * pow(s + t, 2) +
                     2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * t) +
                     pow(mrho, 2) * ((-2 + delta) * s + (-6 + delta) * t +
                                     8 * C4 * pow(s + t, 2)))) +
           eta1 *
               (2 * delta * t * pow(s + t, 2) -
                pow(mrho, 2) *
                    ((-2 + delta) * pow(s, 2) + 2 * (-2 + 3 * delta) * s * t +
                     (-6 + 5 * delta) * pow(t, 2) +
                     8 * C4 * t * pow(s + t, 2)) +
                pow(mrho, 4) *
                    (8 * C4 * t * (2 * s + t) + (-2 + delta) * (s + 3 * t)) +
                4 * pow(mpion, 4) *
                    (6 * C4 * pow(mrho, 4) + 2 * delta * t +
                     pow(mrho, 2) * (1 - 2 * delta - 8 * C4 * t)) -
                2 * pow(mpion, 2) *
                    (4 * delta * t * (s + t) -
                     pow(mrho, 2) * (-2 * s + 3 * delta * s - 6 * t +
                                     7 * delta * t + 16 * C4 * t * (s + t)) +
                     2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * (s + 2 * t))))) /
              (pow(mrho, 2) * (-pow(ma1, 2) + t)))) /
            (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
             2 * pow(mpion, 2) * (pow(mrho, 2) + s)))) /
      (512. * Pi);

  return diff_sigma / spin_deg_factor;
}

// C12
double TwoToTwo_Diff_Xsections::diff_xsection_C12(double t, double s) {
  double diff_sigma;
  double spin_deg_factor = 3.0;

  diff_sigma =
      (pow(Const, 2) * pow(ghat, 4) *
       ((-0.25 * pow(-2 + delta, 2) * pow(mpion, 2) *
         (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
          2 * pow(mpion, 2) * (pow(mrho, 2) + s))) /
            (pow(mrho, 2) * pow(pow(mpion, 2) - s, 2)) -
        (0.0625 * (eta1 - eta2) * (-pow(ma1, 2) + s) *
         (2 * pow(mrho, 2) +
          delta * (-2 * pow(mpion, 2) - pow(mrho, 2) + s + t)) *
         (-(eta2 * (s - t) *
            (4 * pow(mpion, 4) + s * (4 * pow(mrho, 2) + s - t) -
             pow(mpion, 2) * (3 * s + t))) +
          eta1 * (8 * pow(mpion, 6) + pow(s, 3) + pow(s, 2) * t +
                  5 * s * pow(t, 2) + pow(t, 3) + 2 * pow(mrho, 4) * (-s + t) +
                  pow(mrho, 2) * (s - 3 * t) * (s + t) +
                  2 * pow(mpion, 2) * (2 * pow(mrho, 2) - s - t) * (s + t) -
                  2 * pow(mpion, 4) * (2 * pow(mrho, 2) + s + 3 * t)))) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (-2 * pow(mpion, 2) + s + t)) -
        (0.0625 *
         pow(-2. * pow(mrho, 2) +
                 delta * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t),
             2) *
         (8. * pow(mpion, 6) + 4. * pow(mrho, 6) + pow(s, 3) +
          pow(mrho, 4) * (-4. * s - 4. * t) +
          pow(mpion, 4) * (-4. * pow(mrho, 2) - 4. * s - 4. * t) +
          3. * pow(s, 2) * t + 3. * s * pow(t, 2) + pow(t, 3) +
          pow(mrho, 2) * (-3. * pow(s, 2) + 2. * s * t - 3. * pow(t, 2)) +
          pow(mpion, 2) *
              (-8. * pow(mrho, 4) - 2. * pow(s, 2) - 4. * s * t -
               2. * pow(t, 2) + pow(mrho, 2) * (4. * s + 4. * t)))) /
            (pow(mrho, 6) * pow(2. * pow(mpion, 2) - 1. * s - 1. * t, 2)) +
        (0.125 * (-2 + delta) * (eta1 - eta2) * (-pow(ma1, 2) + s) *
         (-(eta2 * (pow(mpion, 2) + s) *
            (-pow(mpion, 4) + pow(mpion, 2) * (pow(mrho, 2) - 2 * s) +
             s * (-pow(mrho, 2) + s + 2 * t))) +
          eta1 * (-4 * pow(mpion, 6) +
                  s * (-pow(mrho, 2) + s) * (-pow(mrho, 2) + s + t) +
                  pow(mpion, 4) * (3 * pow(mrho, 2) + s + t) -
                  pow(mpion, 2) * (pow(mrho, 4) + 2 * s * (s - t) +
                                   pow(mrho, 2) * (-s + t))))) /
            ((pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) *
             (-pow(mpion, 2) + s)) +
        (0.03125 * pow(eta1 - eta2, 2) *
         (pow(eta2, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 4) * (pow(pow(mrho, 2) + 2 * s, 2) - 2 * s * t) +
               pow(s, 2) * (pow(pow(mrho, 2) + s, 2) +
                            2 * (-pow(mrho, 2) + s) * t + 2 * pow(t, 2)) -
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * (2 * s - t) +
                    2 * s * (s + t))) -
          2 * eta1 * eta2 *
              (pow(mpion, 8) -
               pow(mpion, 4) * (pow(mrho, 4) + 2 * pow(mrho, 2) * s +
                                2 * s * (-2 * s + t)) +
               2 * pow(mpion, 2) * s *
                   (pow(mrho, 4) + pow(mrho, 2) * (s + t) - 2 * s * (s + t)) +
               pow(s, 2) * (-pow(mrho, 4) + pow(s, 2) - 2 * pow(mrho, 2) * t +
                            2 * t * (s + t))) +
          pow(eta1, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 4) * (3 * pow(mrho, 4) + 2 * s * (2 * s - t) +
                                2 * pow(mrho, 2) * (-3 * s + t)) -
               2 * pow(mpion, 2) * (-pow(mrho, 2) + s) *
                   (2 * s * (s + t) - pow(mrho, 2) * (2 * s + t)) +
               s * (-pow(mrho, 2) + s) *
                   (pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                    pow(mrho, 2) * (s + 2 * t))))) /
            (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2)) +
        (0.5 *
         (-2. * pow(mrho, 2) +
          delta * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t)) *
         (delta * (1. * pow(mpion, 6) + 0.5 * pow(mrho, 6) +
                   0.0625 * pow(s, 3) + pow(mrho, 4) * (-0.5 * s - 0.5 * t) +
                   pow(mpion, 4) * (-0.5 * pow(mrho, 2) - 0.75 * s - 0.25 * t) +
                   0.3125 * pow(s, 2) * t + 0.4375 * s * pow(t, 2) +
                   0.1875 * pow(t, 3) +
                   pow(mpion, 2) * (-1. * pow(mrho, 4) +
                                    pow(mrho, 2) * (0.375 * s + 0.625 * t) +
                                    (-0.5 * s - 0.5 * t) * t) +
                   pow(mrho, 2) * (-0.3125 * pow(s, 2) + 0.25 * s * t -
                                   0.4375 * pow(t, 2))) +
          pow(mrho, 2) *
              (-0.125 * pow(s, 2) + C4 * pow(mrho, 4) * (1. * s - 1. * t) +
               0.125 * pow(t, 2) +
               pow(mpion, 2) * ((0.25 - 1. * C4 * pow(mrho, 2)) * s +
                                (-0.25 + 1. * C4 * pow(mrho, 2)) * t) +
               pow(mrho, 2) * (-0.5 * s + 0.5 * C4 * pow(s, 2) +
                               t * (0.5 - 0.5 * C4 * t))))) /
            (pow(mrho, 6) * (1. * pow(mpion, 2) - 0.5 * s - 0.5 * t)) +
        (pow(delta, 2) *
             (-0.5 * pow(mpion, 6) - 0.0625 * pow(mrho, 6) +
              pow(mpion, 4) * (1. * pow(mrho, 2) + 0.5 * s) +
              pow(mrho, 4) * (-0.125 * s - 0.125 * t) +
              t * (-0.125 * pow(s, 2) - 0.25 * s * t - 0.125 * pow(t, 2)) +
              pow(mpion, 2) * (1.25 * pow(mrho, 4) - 0.125 * pow(s, 2) +
                               pow(mrho, 2) * (-0.875 * s - 1.125 * t) +
                               0.25 * s * t + 0.375 * pow(t, 2)) +
              pow(mrho, 2) *
                  (0.3125 * pow(s, 2) + 0.25 * s * t + 0.4375 * pow(t, 2))) +
         delta * pow(mrho, 2) *
             (pow(mpion, 4) * (-0.5 - 4. * C4 * pow(mrho, 2)) +
              (-0.25 * s - 0.25 * t) * t +
              pow(mrho, 4) * (-0.75 + 0.5 * C4 * s + 2.5 * C4 * t) +
              pow(mrho, 2) * (-1.5 * C4 * pow(s, 2) + s * (1.25 - 2. * C4 * t) +
                              t * (0.25 - 0.5 * C4 * t)) +
              pow(mpion, 2) *
                  (-3. * C4 * pow(mrho, 4) + 0.25 * s + 0.75 * t +
                   pow(mrho, 2) * (-1.5 + 5. * C4 * s + 3. * C4 * t))) +
         pow(mrho, 6) *
             (0.75 +
              C4 * (8. * C4 * pow(mpion, 4) + 2. * C4 * pow(s, 2) +
                    pow(mpion, 2) * (6. - 8. * C4 * s - 8. * C4 * t) +
                    t * (-3. + 2. * C4 * t) + s * (-3. + 4. * C4 * t)))) /
            pow(mrho, 6) +
        (0.0625 * (eta1 - eta2) * (-pow(ma1, 2) + s) *
         (-(eta2 *
            (2 * pow(mpion, 4) *
                 (-4 * C4 * pow(mrho, 2) * (pow(mrho, 2) + 4 * s) +
                  delta * (pow(mrho, 2) + 6 * s - 2 * t)) +
             pow(mpion, 2) *
                 (2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * s) +
                  delta * (-11 * pow(s, 2) - 6 * s * t + pow(t, 2)) +
                  pow(mrho, 2) * ((-10 + delta) * s - (-2 + delta) * t +
                                  32 * C4 * s * (s + t))) +
             s * (-2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * s) +
                  delta * (3 * pow(s, 2) + 2 * s * t + 3 * pow(t, 2)) +
                  pow(mrho, 2) * (3 * (2 + delta) * s + (2 - 5 * delta) * t -
                                  8 * C4 * pow(s + t, 2))))) +
          eta1 *
              (8 * delta * pow(mpion, 6) +
               2 * pow(mpion, 4) *
                   (12 * C4 * pow(mrho, 4) -
                    2 * pow(mrho, 2) * (-1 + 3 * delta + 8 * C4 * s) +
                    3 * delta * (s - t)) +
               delta * (3 * pow(s, 3) + 5 * pow(s, 2) * t + 7 * s * pow(t, 2) +
                        pow(t, 3)) -
               2 * pow(mrho, 2) *
                   ((-3 + 2 * delta) * pow(s, 2) +
                    2 * (-1 + 2 * delta) * s * t +
                    (-1 + 2 * delta) * pow(t, 2) + 4 * C4 * s * pow(s + t, 2)) +
               pow(mrho, 4) * ((-6 + delta) * s + (-2 + 3 * delta) * t +
                               8 * C4 * s * (s + 2 * t)) -
               2 * pow(mpion, 2) *
                   (delta * (s + t) * (5 * s + t) -
                    pow(mrho, 2) * (-6 * s + 9 * delta * s - 2 * t +
                                    5 * delta * t + 16 * C4 * s * (s + t)) +
                    2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * (2 * s + t)))))) /
            (pow(mrho, 2) *
             (pow(Gammaa1, 2) * pow(ma1, 2) + pow(pow(ma1, 2) - s, 2))) +
        (2 *
         ((0.0625 * (-2. + delta) *
           (-2. * pow(mrho, 2) +
            delta * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t)) *
           (2. * pow(mpion, 6) + 1. * pow(mrho, 6) +
            pow(mpion, 4) * (-3. * pow(mrho, 2) - 2. * s) +
            pow(mrho, 4) * (-1.5 * s - 1.5 * t) +
            pow(mrho, 2) * (0.5 * s + 0.5 * t) * t +
            s * (0.5 * pow(s, 2) + 1. * s * t + 0.5 * pow(t, 2)) +
            pow(mpion, 2) *
                (-1. * pow(mrho, 4) - 0.5 * pow(s, 2) - 1. * s * t -
                 0.5 * pow(t, 2) + pow(mrho, 2) * (-0.5 * s + 2.5 * t)))) /
              ((pow(mpion, 2) - 1. * s) *
               (1. * pow(mpion, 2) - 0.5 * s - 0.5 * t)) +
          (0.0625 * (-2 + delta) *
           (delta *
                (6 * pow(mpion, 6) -
                 pow(mpion, 4) * (9 * (pow(mrho, 2) + s) + t) -
                 pow(mpion, 2) * (2 * pow(mrho, 4) - 3 * pow(s, 2) +
                                  pow(mrho, 2) * (5 * s - 7 * t) + pow(t, 2)) +
                 (pow(mrho, 2) - s - t) *
                     (3 * pow(mrho, 4) - s * t - pow(mrho, 2) * (3 * s + t))) +
            2 * pow(mrho, 2) *
                (pow(mpion, 4) * (1 + 4 * C4 * pow(mrho, 2)) +
                 pow(mrho, 4) * (-1 + 4 * C4 * s) + s * t -
                 pow(mpion, 2) * (4 * C4 * pow(mrho, 4) + s -
                                  2 * pow(mrho, 2) * (1 + 4 * C4 * s) + t) +
                 pow(mrho, 2) * (t + s * (1 - 4 * C4 * (s + 2 * t)))))) /
              (-pow(mpion, 2) + s))) /
            pow(mrho, 4))) /
      (16. * Pi *
       (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
        2 * pow(mpion, 2) * (pow(mrho, 2) + s)));

  return diff_sigma / spin_deg_factor;
}

// C13
double TwoToTwo_Diff_Xsections::diff_xsection_C13(double t, double s) {
  double diff_sigma;
  double spin_deg_factor = 3.0;

  diff_sigma =
      (pow(Const, 2) * pow(ghat, 4) *
       ((-0.25 * pow(-2 + delta, 2) * pow(mpion, 2) *
         (pow(mpion, 4) + pow(pow(mrho, 2) - t, 2) -
          2 * pow(mpion, 2) * (pow(mrho, 2) + t))) /
            (pow(mrho, 2) * pow(pow(mpion, 2) - t, 2)) -
        (0.0625 * (eta1 - eta2) *
         (2 * pow(mrho, 2) +
          delta * (-2 * pow(mpion, 2) - pow(mrho, 2) + s + t)) *
         (eta1 * (8 * pow(mpion, 6) + pow(s, 3) + 2 * pow(mrho, 4) * (s - t) +
                  5 * pow(s, 2) * t + s * pow(t, 2) + pow(t, 3) +
                  2 * pow(mpion, 2) * (2 * pow(mrho, 2) - s - t) * (s + t) -
                  pow(mrho, 2) * (3 * s - t) * (s + t) -
                  2 * pow(mpion, 4) * (2 * pow(mrho, 2) + 3 * s + t)) +
          eta2 * (s - t) *
              (4 * pow(mpion, 4) + t * (4 * pow(mrho, 2) - s + t) -
               pow(mpion, 2) * (s + 3 * t)))) /
            (pow(mrho, 2) * (-pow(ma1, 2) + t) * (-2 * pow(mpion, 2) + s + t)) -
        (0.0625 *
         pow(-2. * pow(mrho, 2) +
                 delta * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t),
             2) *
         (8. * pow(mpion, 6) + 4. * pow(mrho, 6) + pow(s, 3) +
          pow(mrho, 4) * (-4. * s - 4. * t) +
          pow(mpion, 4) * (-4. * pow(mrho, 2) - 4. * s - 4. * t) +
          3. * pow(s, 2) * t + 3. * s * pow(t, 2) + pow(t, 3) +
          pow(mrho, 2) * (-3. * pow(s, 2) + 2. * s * t - 3. * pow(t, 2)) +
          pow(mpion, 2) *
              (-8. * pow(mrho, 4) - 2. * pow(s, 2) - 4. * s * t -
               2. * pow(t, 2) + pow(mrho, 2) * (4. * s + 4. * t)))) /
            (pow(mrho, 6) * pow(2. * pow(mpion, 2) - 1. * s - 1. * t, 2)) +
        (0.125 * (-2 + delta) * (eta1 - eta2) *
         (eta2 * (pow(mpion, 2) + t) *
              (pow(mpion, 4) - pow(mpion, 2) * (pow(mrho, 2) - 2 * t) +
               (pow(mrho, 2) - 2 * s - t) * t) +
          eta1 * (-4 * pow(mpion, 6) +
                  (pow(mrho, 2) - t) * (pow(mrho, 2) - s - t) * t +
                  pow(mpion, 4) * (3 * pow(mrho, 2) + s + t) -
                  pow(mpion, 2) * (pow(mrho, 4) + pow(mrho, 2) * (s - t) +
                                   2 * t * (-s + t))))) /
            ((-pow(ma1, 2) + t) * (-pow(mpion, 2) + t)) +
        (0.03125 * pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(mpion, 8) -
               pow(mpion, 4) *
                   (pow(mrho, 4) + 2 * (pow(mrho, 2) + s) * t - 4 * pow(t, 2)) +
               pow(t, 2) * (-pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                            2 * pow(s, 2) + 2 * s * t + pow(t, 2)) +
               2 * pow(mpion, 2) * t *
                   (pow(mrho, 4) + pow(mrho, 2) * (s + t) - 2 * t * (s + t))) +
          pow(eta2, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 4) *
                   (pow(mrho, 4) + 4 * pow(mrho, 2) * t - 2 * (s - 2 * t) * t) +
               pow(t, 2) * (pow(mrho, 4) + 2 * pow(s, 2) + 2 * s * t +
                            pow(t, 2) + 2 * pow(mrho, 2) * (-s + t)) -
               2 * pow(mpion, 2) * t *
                   (pow(mrho, 4) - pow(mrho, 2) * (s - 2 * t) +
                    2 * t * (s + t))) +
          pow(eta1, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
               pow(mpion, 4) *
                   (3 * pow(mrho, 4) + 2 * pow(mrho, 2) * (s - 3 * t) -
                    2 * (s - 2 * t) * t) +
               t * (-pow(mrho, 2) + t) *
                   (2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                    pow(mrho, 2) * (2 * s + t)) -
               2 * pow(mpion, 2) * (-pow(mrho, 2) + t) *
                   (2 * t * (s + t) - pow(mrho, 2) * (s + 2 * t))))) /
            pow(pow(ma1, 2) - t, 2) -
        (0.5 *
         (-2. * pow(mrho, 2) +
          delta * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t)) *
         (delta *
              (-1. * pow(mpion, 6) - 0.5 * pow(mrho, 6) - 0.1875 * pow(s, 3) +
               pow(mpion, 2) * (1. * pow(mrho, 4) +
                                pow(mrho, 2) * (-0.625 * s - 0.375 * t) +
                                s * (0.5 * s + 0.5 * t)) +
               pow(mrho, 4) * (0.5 * s + 0.5 * t) +
               pow(mpion, 4) * (0.5 * pow(mrho, 2) + 0.25 * s + 0.75 * t) -
               0.4375 * pow(s, 2) * t - 0.3125 * s * pow(t, 2) -
               0.0625 * pow(t, 3) +
               pow(mrho, 2) *
                   (0.4375 * pow(s, 2) - 0.25 * s * t + 0.3125 * pow(t, 2))) +
          pow(mrho, 2) *
              (-0.125 * pow(s, 2) + C4 * pow(mrho, 4) * (1. * s - 1. * t) +
               0.125 * pow(t, 2) +
               pow(mpion, 2) * ((0.25 - 1. * C4 * pow(mrho, 2)) * s +
                                (-0.25 + 1. * C4 * pow(mrho, 2)) * t) +
               pow(mrho, 2) * (-0.5 * s + 0.5 * C4 * pow(s, 2) +
                               t * (0.5 - 0.5 * C4 * t))))) /
            (pow(mrho, 6) * (1. * pow(mpion, 2) - 0.5 * s - 0.5 * t)) +
        (pow(delta, 2) *
             (-0.5 * pow(mpion, 6) - 0.0625 * pow(mrho, 6) +
              pow(mrho, 4) * (-0.125 * s - 0.125 * t) +
              pow(mpion, 4) * (1. * pow(mrho, 2) + 0.5 * t) +
              s * (-0.125 * pow(s, 2) - 0.25 * s * t - 0.125 * pow(t, 2)) +
              pow(mpion, 2) * (1.25 * pow(mrho, 4) + 0.375 * pow(s, 2) +
                               pow(mrho, 2) * (-1.125 * s - 0.875 * t) +
                               0.25 * s * t - 0.125 * pow(t, 2)) +
              pow(mrho, 2) *
                  (0.4375 * pow(s, 2) + 0.25 * s * t + 0.3125 * pow(t, 2))) +
         pow(mrho, 6) *
             (0.75 + C4 * (8. * C4 * pow(mpion, 4) + 2. * C4 * pow(s, 2) +
                           pow(mpion, 2) * (6. - 8. * C4 * s - 8. * C4 * t) +
                           t * (-3. + 2. * C4 * t) + s * (-3. + 4. * C4 * t))) +
         delta * pow(mrho, 2) *
             (pow(mpion, 4) * (-0.5 - 4. * C4 * pow(mrho, 2)) +
              s * (-0.25 * s - 0.25 * t) +
              pow(mrho, 4) * (-0.75 + 2.5 * C4 * s + 0.5 * C4 * t) +
              pow(mrho, 2) * (-0.5 * C4 * pow(s, 2) + s * (0.25 - 2. * C4 * t) +
                              t * (1.25 - 1.5 * C4 * t)) +
              pow(mpion, 2) *
                  (-3. * C4 * pow(mrho, 4) + 0.75 * s + 0.25 * t +
                   pow(mrho, 2) * (-1.5 + 3. * C4 * s + 5. * C4 * t)))) /
            pow(mrho, 6) +
        (2 * ((0.0625 * (-2. + delta) *
               (-2. * pow(mrho, 2) +
                delta * (2. * pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t)) *
               (2. * pow(mpion, 6) + 1. * pow(mrho, 6) +
                pow(mpion, 4) * (-3. * pow(mrho, 2) - 2. * t) +
                pow(mrho, 4) * (-1.5 * s - 1.5 * t) +
                pow(mrho, 2) * s * (0.5 * s + 0.5 * t) +
                pow(mpion, 2) * (-1. * pow(mrho, 4) - 0.5 * pow(s, 2) +
                                 pow(mrho, 2) * (2.5 * s - 0.5 * t) -
                                 1. * s * t - 0.5 * pow(t, 2)) +
                t * (0.5 * pow(s, 2) + 1. * s * t + 0.5 * pow(t, 2)))) /
                  ((pow(mpion, 2) - 1. * t) *
                   (1. * pow(mpion, 2) - 0.5 * s - 0.5 * t)) +
              (0.0625 * (-2 + delta) *
               (6 * delta * pow(mpion, 6) + delta * s * t * (s + t) +
                pow(mrho, 6) * (-2 + 3 * delta + 8 * C4 * t) -
                pow(mpion, 4) * ((-2 + 9 * delta) * pow(mrho, 2) -
                                 8 * C4 * pow(mrho, 4) + delta * (s + 9 * t)) -
                2 * pow(mrho, 4) *
                    (t * (-1 + 3 * delta + 4 * C4 * t) +
                     s * (-1 + 2 * delta + 8 * C4 * t)) -
                pow(mpion, 2) * (8 * C4 * pow(mrho, 6) +
                                 2 * pow(mrho, 4) * (-2 + delta - 8 * C4 * t) +
                                 pow(mrho, 2) * ((2 - 7 * delta) * s +
                                                 (2 + 5 * delta) * t) +
                                 delta * (pow(s, 2) - 3 * pow(t, 2))) +
                pow(mrho, 2) * (2 * s * t + delta * (pow(s, 2) + 3 * s * t +
                                                     3 * pow(t, 2))))) /
                  (-pow(mpion, 2) + t))) /
            pow(mrho, 4) +
        (0.0625 * (eta1 - eta2) *
         (-(eta2 *
            (-2 * pow(mpion, 4) *
                 (4 * C4 * pow(mrho, 2) * (pow(mrho, 2) + 4 * t) -
                  delta * (pow(mrho, 2) - 2 * s + 6 * t)) +
             pow(mpion, 2) *
                 (2 * pow(mrho, 4) * (-2 + delta + 8 * C4 * t) +
                  delta * (pow(s, 2) - 6 * s * t - 11 * pow(t, 2)) +
                  pow(mrho, 2) * (-((-2 + delta) * s) + (-10 + delta) * t +
                                  32 * C4 * t * (s + t))) +
             t * (-2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * t) +
                  delta * (3 * pow(s, 2) + 2 * s * t + 3 * pow(t, 2)) +
                  pow(mrho, 2) * ((2 - 5 * delta) * s + 3 * (2 + delta) * t -
                                  8 * C4 * pow(s + t, 2))))) +
          eta1 *
              (8 * delta * pow(mpion, 6) +
               delta * (pow(s, 3) + 7 * pow(s, 2) * t + 5 * s * pow(t, 2) +
                        3 * pow(t, 3)) -
               2 * pow(mrho, 2) *
                   ((-1 + 2 * delta) * pow(s, 2) +
                    2 * (-1 + 2 * delta) * s * t +
                    (-3 + 2 * delta) * pow(t, 2) + 4 * C4 * t * pow(s + t, 2)) +
               pow(mpion, 4) *
                   (24 * C4 * pow(mrho, 4) + 6 * delta * (-s + t) -
                    4 * pow(mrho, 2) * (-1 + 3 * delta + 8 * C4 * t)) +
               pow(mrho, 4) * (t * (-6 + delta + 8 * C4 * t) +
                               s * (-2 + 3 * delta + 16 * C4 * t)) -
               2 * pow(mpion, 2) *
                   (delta * (s + t) * (s + 5 * t) -
                    pow(mrho, 2) * (-2 * s + 5 * delta * s - 6 * t +
                                    9 * delta * t + 16 * C4 * t * (s + t)) +
                    2 * pow(mrho, 4) * (-2 + delta + 4 * C4 * (s + 2 * t)))))) /
            (pow(mrho, 2) * (-pow(ma1, 2) + t)))) /
      (16. * Pi *
       (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
        2 * pow(mpion, 2) * (pow(mrho, 2) + s)));

  return diff_sigma / spin_deg_factor;
}

// C14
double TwoToTwo_Diff_Xsections::diff_xsection_C14(double t, double s) {
  double diff_sigma;
  double spin_deg_factor = 3.0;

  diff_sigma =
      (pow(Const, 2) * pow(g_POR, 4) *
       ((0.125 *
         (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
          pow(mpion, 4) * (pow(mrho, 4) - 2 * (s - 2 * t) * t) +
          pow(t, 2) * (pow(mrho, 4) + 2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                       2 * pow(mrho, 2) * (s + t)) -
          2 * pow(mpion, 2) * t *
              (pow(mrho, 4) + 2 * t * (s + t) - pow(mrho, 2) * (s + 2 * t)))) /
            pow(pow(momega, 2) - t, 2) +
        (((0.25 * (-pow(momega, 2) + s) *
           (pow(mpion, 8) -
            4 * pow(mpion, 2) * s * t * (-pow(mrho, 2) + s + t) +
            pow(mpion, 4) * (2 * s * t - pow(mrho, 2) * (s + t)) +
            s * t *
                (pow(s, 2) + 3 * s * t + pow(t, 2) - pow(mrho, 2) * (s + t)))) /
              (-pow(momega, 2) + t) +
          0.125 * (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
                   pow(mpion, 4) * (pow(mrho, 4) + 4 * pow(s, 2) - 2 * s * t) +
                   pow(s, 2) * (pow(mrho, 4) + pow(s, 2) + 2 * s * t +
                                2 * pow(t, 2) - 2 * pow(mrho, 2) * (s + t)) -
                   2 * pow(mpion, 2) * s *
                       (pow(mrho, 4) + 2 * s * (s + t) -
                        pow(mrho, 2) * (2 * s + t)))) *
         HeavisideTheta(-momega + sqrt(s))) /
            pow(pow(momega, 2) - s, 2))) /
      (16. * Pi *
       (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
        2 * pow(mpion, 2) * (pow(mrho, 2) + s)));

  return diff_sigma / spin_deg_factor;
}

// C15
double TwoToTwo_Diff_Xsections::diff_xsection_C15(double t, double s) {
  double diff_sigma;
  double spin_deg_factor = 3.0;

  diff_sigma =
      HeavisideTheta(-momega + sqrt(s)) *
      ((0.0024867959858108648 * pow(Const, 2) * pow(g_POR, 4) *
        (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
         pow(mpion, 4) * (pow(mrho, 4) + 4 * pow(s, 2) - 2 * s * t) +
         pow(s, 2) * (pow(mrho, 4) + pow(s, 2) + 2 * s * t + 2 * pow(t, 2) -
                      2 * pow(mrho, 2) * (s + t)) -
         2 * pow(mpion, 2) * s *
             (pow(mrho, 4) + 2 * s * (s + t) - pow(mrho, 2) * (2 * s + t)))) /
       ((pow(pow(momega, 2) - s, 2) *
         (pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
          2 * pow(mpion, 2) * (pow(mrho, 2) + s)))));

  return diff_sigma / spin_deg_factor;
}

// C16
double TwoToTwo_Diff_Xsections::diff_xsection_C16(double t, double s) {
  double diff_sigma;
  double spin_deg_factor = 3.0;

  diff_sigma =
      (0.0024867959858108648 * pow(Const, 2) * pow(g_POR, 4) *
       (pow(mpion, 8) - 2 * pow(mpion, 6) * pow(mrho, 2) +
        pow(mpion, 4) * (pow(mrho, 4) - 2 * (s - 2 * t) * t) +
        pow(t, 2) * (pow(mrho, 4) + 2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                     2 * pow(mrho, 2) * (s + t)) -
        2 * pow(mpion, 2) * t *
            (pow(mrho, 4) + 2 * t * (s + t) - pow(mrho, 2) * (s + 2 * t)))) /
      ((pow(mpion, 4) + pow(pow(mrho, 2) - s, 2) -
        2 * pow(mpion, 2) * (pow(mrho, 2) + s)) *
       pow(pow(momega, 2) - t, 2));

  return diff_sigma / spin_deg_factor;
}

// C21
double TwoToTwo_Diff_Xsections::diff_xsection_C21(double t, double s) {
  double diff_sigma;
  double spin_deg_factor = 1.0;

  diff_sigma =
      ((pow(Const, 2) * pow(ghat, 4) *
        ((0.25 *
          (32 * pow(C4, 2) * pow(mrho, 8) + 2 * pow(delta, 2) * pow(s, 2) +
           8 * C4 * pow(mrho, 6) * (-6 + delta - 8 * C4 * s) +
           2 * delta * pow(mrho, 2) * s * (-6 + delta - 8 * C4 * s) +
           pow(mrho, 4) * (12 - pow(delta, 2) + 8 * C4 * (6 + delta) * s +
                           32 * pow(C4, 2) * pow(s, 2)))) /
             pow(mrho, 4) -
         (0.25 * pow(-2 + delta, 2) * pow(mpion, 2) *
          (pow(mpion, 4) + pow(pow(mrho, 2) - t, 2) -
           2 * pow(mpion, 2) * (pow(mrho, 2) + t))) /
             (pow(mrho, 2) * pow(pow(mpion, 2) - t, 2)) -
         (0.25 * pow(-2 + delta, 2) * pow(mpion, 2) *
          (pow(mpion, 4) + pow(s + t, 2) -
           2 * pow(mpion, 2) * (2 * pow(mrho, 2) + s + t))) /
             (pow(mrho, 2) * pow(pow(mpion, 2) + pow(mrho, 2) - s - t, 2)) +
         (0.125 * (-2 + delta) * (eta1 - eta2) *
          (pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t) *
          (eta1 * (2 * pow(mpion, 2) - s) +
           eta2 * (-3 * pow(mpion, 2) - pow(mrho, 2) + s + t)) *
          (pow(mpion, 4) + t * (-pow(mrho, 2) + 2 * s + t) -
           pow(mpion, 2) * (pow(mrho, 2) + 2 * t))) /
             ((-pow(mpion, 2) + t) *
              (pow(Gammaa1, 2) * pow(ma1, 2) +
               pow(pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t,
                   2))) -
         (0.125 * (-2 + delta) * (eta1 - eta2) *
          (eta1 * (-2 * pow(mpion, 2) + s) + eta2 * (pow(mpion, 2) + t)) *
          (-pow(mpion, 4) + pow(s, 2) - pow(t, 2) + pow(mrho, 2) * (-s + t) +
           pow(mpion, 2) * (pow(mrho, 2) - 2 * s + 2 * t))) /
             ((pow(mpion, 2) + pow(mrho, 2) - s - t) * (-pow(ma1, 2) + t)) +
         (0.25 * (-2. + delta) *
          (pow(mpion, 4) * (2. + delta - 8. * C4 * pow(mrho, 2)) +
           8. * C4 * pow(mrho, 4) * t +
           t * ((2. + 3. * delta) * s + (2. + delta) * t) +
           pow(mrho, 2) * (s * (2. - 1. * delta - 16. * C4 * t) +
                           t * (-2. - 1. * delta - 8. * C4 * t)) +
           pow(mpion, 2) * (8. * C4 * pow(mrho, 4) + (-2. + delta) * s +
                            (-4. - 2. * delta) * t +
                            pow(mrho, 2) * (-6. + delta + 16. * C4 * t)))) /
             (pow(mrho, 2) * (pow(mpion, 2) - 1. * t)) -
         (0.125 * (-2 + delta) * (eta1 - eta2) *
          (pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t) *
          (-(eta2 * (3 * pow(mpion, 2) + pow(mrho, 2) - s - t) *
             (pow(mpion, 4) + (pow(mrho, 2) - s - t) * (s - t) -
              pow(mpion, 2) * (pow(mrho, 2) - 2 * s + 2 * t))) +
           eta1 *
               (2 * pow(mpion, 6) +
                pow(mpion, 4) * (-2 * pow(mrho, 2) + 5 * s - 4 * t) +
                s * (s + t) * (-pow(mrho, 2) + s + t) +
                pow(mpion, 2) * (2 * pow(mrho, 4) + pow(mrho, 2) * (s - 2 * t) -
                                 2 * (2 * s - t) * (s + t))))) /
             ((-pow(mpion, 2) - pow(mrho, 2) + s + t) *
              (pow(Gammaa1, 2) * pow(ma1, 2) +
               pow(pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t,
                   2))) +
         (0.25 * (-2. + delta) * (eta1 - 1. * eta2) *
          (eta2 * (-0.5 * pow(mpion, 6) +
                   pow(mpion, 4) * (0.5 * pow(mrho, 2) + 0.5 * t) +
                   pow(mpion, 2) * (1. * pow(mrho, 2) - 1. * s + 0.5 * t) * t +
                   (0.5 * pow(mrho, 2) - 1. * s - 0.5 * t) * pow(t, 2)) +
           eta1 * (1. * pow(mpion, 6) +
                   pow(mpion, 4) * (-1. * pow(mrho, 2) + 0.5 * s - 2. * t) +
                   s * (-0.5 * pow(mrho, 2) + 0.5 * t) * t +
                   pow(mpion, 2) *
                       (1. * pow(mrho, 4) + pow(mrho, 2) * (-0.5 * s - 1. * t) +
                        t * (1. * s + 1. * t))))) /
             ((pow(ma1, 2) - 1. * t) * (-1. * pow(mpion, 2) + t)) +
         (0.03125 * pow(eta1 - eta2, 2) *
          (-2 * eta1 * eta2 *
               (pow(mpion, 8) - 4 * pow(mpion, 6) * t +
                pow(t, 2) * (-pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                             2 * pow(s, 2) + 2 * s * t + pow(t, 2)) -
                2 * pow(mpion, 2) * t *
                    (-2 * pow(mrho, 4) + pow(mrho, 2) * s + 2 * t * (s + t)) +
                pow(mpion, 4) * (-pow(mrho, 4) + 2 * t * (s + 3 * t))) +
           pow(eta2, 2) *
               (pow(mpion, 8) - 2 * pow(mpion, 6) * (pow(mrho, 2) + 2 * t) +
                pow(t, 2) * (pow(mrho, 4) + 2 * pow(s, 2) + 2 * s * t +
                             pow(t, 2) + 2 * pow(mrho, 2) * (-s + t)) -
                2 * pow(mpion, 2) * t *
                    (2 * t * (s + t) + pow(mrho, 2) * (s + 3 * t)) +
                pow(mpion, 4) * (pow(mrho, 4) + 6 * pow(mrho, 2) * t +
                                 2 * t * (s + 3 * t))) +
           pow(eta1, 2) *
               (pow(mpion, 8) + 2 * pow(mpion, 6) * (pow(mrho, 2) - 2 * t) -
                2 * pow(mpion, 2) * (pow(mrho, 2) - s - t) *
                    (pow(mrho, 4) + pow(mrho, 2) * t - 2 * pow(t, 2)) +
                t * (-pow(mrho, 2) + t) *
                    (2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                     pow(mrho, 2) * (2 * s + t)) +
                pow(mpion, 4) * (pow(mrho, 4) - 2 * pow(mrho, 2) * (s + 3 * t) +
                                 2 * t * (s + 3 * t))))) /
             pow(pow(ma1, 2) - t, 2) +
         ((0.25 * pow(-2 + delta, 2) * (2 * pow(mpion, 2) - s) *
           (pow(mpion, 4) + pow(mrho, 2) * (s - t) + t * (s + t) -
            pow(mpion, 2) * (3 * pow(mrho, 2) + s + 2 * t))) /
              ((pow(mpion, 2) - t) * (pow(mpion, 2) + pow(mrho, 2) - s - t)) -
          (0.25 * (-2. + delta) *
           (pow(mpion, 4) * (2. + delta - 8. * C4 * pow(mrho, 2)) -
            2. * delta * pow(s, 2) + 2. * s * t - 1. * delta * s * t +
            2. * pow(t, 2) + delta * pow(t, 2) +
            C4 * pow(mrho, 4) * (-8. * s + 8. * t) +
            pow(mrho, 2) * ((2. + delta) * s + 8. * C4 * pow(s, 2) +
                            t * (-2. - 1. * delta - 8. * C4 * t)) +
            pow(mpion, 2) *
                (8. * C4 * pow(mrho, 4) - 2. * s + 5. * delta * s - 4. * t -
                 2. * delta * t +
                 pow(mrho, 2) * (-6. + delta - 16. * C4 * s + 16. * C4 * t)))) /
              (pow(mpion, 2) + pow(mrho, 2) - 1. * s - 1. * t)) /
             pow(mrho, 2) +
         (0.03125 * pow(eta1 - eta2, 2) *
          (-2 * eta1 * eta2 *
               (pow(mpion, 8) + 4 * pow(mpion, 6) * (pow(mrho, 2) - t) +
                pow(-pow(mrho, 2) + s + t, 2) *
                    (pow(s, 2) + pow(t, 2) - 2 * pow(mrho, 2) * (s + t)) +
                pow(mpion, 4) *
                    (9 * pow(mrho, 4) + 4 * pow(s, 2) + 2 * s * t +
                     6 * pow(t, 2) - 2 * pow(mrho, 2) * (7 * s + 6 * t)) +
                2 * pow(mpion, 2) * (pow(mrho, 2) - s - t) *
                    (2 * pow(mrho, 4) - pow(mrho, 2) * (5 * s + 4 * t) +
                     2 * (pow(s, 2) + pow(t, 2)))) +
           pow(eta2, 2) *
               (pow(mpion, 8) + pow(mpion, 6) * (6 * pow(mrho, 2) - 4 * t) +
                pow(-pow(mrho, 2) + s + t, 2) *
                    (4 * pow(mrho, 4) + pow(s, 2) + pow(t, 2) -
                     4 * pow(mrho, 2) * (s + t)) +
                pow(mpion, 4) *
                    (17 * pow(mrho, 4) + 4 * pow(s, 2) + 2 * s * t +
                     6 * pow(t, 2) - 2 * pow(mrho, 2) * (10 * s + 9 * t)) +
                2 * pow(mpion, 2) * (pow(mrho, 2) - s - t) *
                    (7 * pow(mrho, 4) - pow(mrho, 2) * (8 * s + 7 * t) +
                     2 * (pow(s, 2) + pow(t, 2)))) +
           pow(eta1, 2) *
               (pow(mpion, 8) + 2 * pow(mpion, 6) * (pow(mrho, 2) - 2 * t) +
                (s + t) * (-pow(mrho, 2) + s + t) *
                    (pow(s, 2) + pow(t, 2) - pow(mrho, 2) * (s + t)) +
                pow(mpion, 4) *
                    (5 * pow(mrho, 4) + 4 * pow(s, 2) + 2 * s * t +
                     6 * pow(t, 2) - 2 * pow(mrho, 2) * (5 * s + 3 * t)) -
                2 * pow(mpion, 2) *
                    (2 * pow(mrho, 4) * (s + t) +
                     2 * (s + t) * (pow(s, 2) + pow(t, 2)) -
                     pow(mrho, 2) *
                         (4 * pow(s, 2) + 5 * s * t + 3 * pow(t, 2)))))) /
             (pow(Gammaa1, 2) * pow(ma1, 2) +
              pow(pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t, 2)) +
         (0.0625 * pow(eta1 - eta2, 2) *
          (pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t) *
          (-(pow(eta2, 2) *
             (pow(mpion, 8) + 2 * pow(mpion, 6) * (pow(mrho, 2) - 2 * t) +
              2 * pow(mpion, 2) * t *
                  (pow(pow(mrho, 2) - s, 2) + (3 * pow(mrho, 2) - 2 * s) * t -
                   2 * pow(t, 2)) +
              (pow(mrho, 2) - s - t) * t *
                  (2 * pow(mrho, 4) + pow(s, 2) - s * t - pow(t, 2) +
                   pow(mrho, 2) * (-3 * s + t)) +
              pow(mpion, 4) * (pow(mrho, 4) + 2 * t * (s + 3 * t) -
                               pow(mrho, 2) * (s + 6 * t)))) -
           pow(eta1, 2) *
               (pow(mpion, 8) + 2 * pow(mpion, 6) * (pow(mrho, 2) - 2 * t) +
                (pow(mrho, 2) - s - t) * t *
                    (pow(s, 2) - s * t - pow(t, 2) + pow(mrho, 2) * (s + t)) +
                pow(mpion, 4) * (3 * pow(mrho, 4) + 2 * t * (s + 3 * t) -
                                 pow(mrho, 2) * (5 * s + 6 * t)) +
                2 * pow(mpion, 2) *
                    (-(pow(mrho, 4) * (s + t)) +
                     t * (pow(s, 2) - 2 * s * t - 2 * pow(t, 2)) +
                     pow(mrho, 2) * (pow(s, 2) + 2 * s * t + 3 * pow(t, 2)))) +
           2 * eta1 * eta2 *
               (pow(mpion, 8) + 2 * pow(mpion, 6) * (pow(mrho, 2) - 2 * t) -
                (pow(mrho, 2) - s - t) * t *
                    (pow(mrho, 4) - pow(s, 2) - pow(mrho, 2) * t +
                     t * (s + t)) +
                2 * pow(mpion, 4) *
                    (2 * pow(mrho, 4) + t * (s + 3 * t) -
                     pow(mrho, 2) * (2 * s + 3 * t)) +
                pow(mpion, 2) *
                    (pow(mrho, 6) - 2 * pow(mrho, 4) * (s + 2 * t) +
                     2 * t * (pow(s, 2) - 2 * s * t - 2 * pow(t, 2)) +
                     pow(mrho, 2) *
                         (pow(s, 2) + 2 * s * t + 6 * pow(t, 2)))))) /
             ((-pow(ma1, 2) + t) *
              (pow(Gammaa1, 2) * pow(ma1, 2) +
               pow(pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t,
                   2))) +
         (0.125 * (eta1 - eta2) *
          (eta2 *
               (8 * C4 * pow(mrho, 6) * t - 2 * delta * pow(s, 2) * t +
                pow(mrho, 2) * (-4 * pow(mpion, 4) +
                                (s * (2 + 3 * delta + 8 * C4 * s) - 4 * t) * t +
                                pow(mpion, 2) * (-((-2 + delta) * s) + 8 * t)) +
                pow(mrho, 4) * (8 * C4 * pow(mpion, 4) -
                                pow(mpion, 2) * (-2 + delta + 16 * C4 * t) +
                                t * (-6 + delta + 8 * C4 * (-2 * s + t)))) +
           eta1 *
               (2 * delta * pow(s, 2) * t +
                8 * C4 * pow(mrho, 6) * (-2 * pow(mpion, 2) + t) -
                pow(mrho, 2) * (-4 * pow(mpion, 4) - 4 * pow(t, 2) +
                                2 * pow(mpion, 2) * ((2 + delta) * s + 4 * t) +
                                pow(s, 2) * (-2 + delta + 8 * C4 * t)) +
                pow(mrho, 4) * (-8 * C4 * pow(mpion, 4) + (-2 + delta) * s -
                                4 * t * (1 + 2 * C4 * t) +
                                8 * pow(mpion, 2) * (1 + 2 * C4 * (s + t)))))) /
             (pow(mrho, 2) * (-pow(ma1, 2) + t)) -
         (0.125 * (eta1 - eta2) *
          (pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t) *
          (eta1 *
               (pow(mpion, 4) * (4 * pow(mrho, 2) - 8 * C4 * pow(mrho, 4)) +
                8 * C4 * pow(mrho, 6) * (s + t) -
                2 * delta * pow(s, 2) * (s + t) +
                pow(mrho, 2) * ((6 + delta) * pow(s, 2) + 8 * s * t +
                                4 * pow(t, 2) + 8 * C4 * pow(s, 2) * (s + t)) -
                pow(mrho, 4) *
                    (-((-6 + delta) * s) + 4 * t +
                     8 * C4 * (2 * pow(s, 2) + 2 * s * t + pow(t, 2))) +
                2 * pow(mpion, 2) *
                    (-8 * C4 * pow(mrho, 6) + 2 * delta * pow(s, 2) -
                     pow(mrho, 2) * (s * (6 + delta + 8 * C4 * s) + 4 * t) +
                     4 * pow(mrho, 4) * (1 + 2 * C4 * (2 * s + t)))) +
           eta2 *
               (pow(mpion, 4) * (-4 * pow(mrho, 2) + 8 * C4 * pow(mrho, 4)) -
                (-pow(mrho, 2) + s + t) *
                    (16 * C4 * pow(mrho, 6) - 2 * delta * pow(s, 2) +
                     pow(mrho, 2) * (s * (6 + 3 * delta + 8 * C4 * s) + 4 * t) +
                     pow(mrho, 4) * (-10 + delta - 8 * C4 * (3 * s + t))) +
                pow(mpion, 2) *
                    (32 * C4 * pow(mrho, 6) - 4 * delta * pow(s, 2) +
                     pow(mrho, 2) *
                         (s * (14 + 5 * delta + 16 * C4 * s) + 8 * t) +
                     pow(mrho, 4) *
                         (delta - 2 * (9 + 8 * C4 * (3 * s + t))))))) /
             (pow(mrho, 2) *
              (pow(Gammaa1, 2) * pow(ma1, 2) +
               pow(pow(ma1, 2) - 2 * pow(mpion, 2) - pow(mrho, 2) + s + t,
                   2))))) /
       (16. * Pi * s * (-4 * pow(mpion, 2) + s)));

  return diff_sigma / spin_deg_factor;
}

// C22
double TwoToTwo_Diff_Xsections::diff_xsection_C22(double t, double s) {
  double diff_sigma;
  double spin_deg_factor = 1.0;

  diff_sigma =
      (pow(Const, 2) * pow(ghat, 4) *
       (0.75 +
        C4 * (2. * C4 * pow(mrho, 4) + pow(mrho, 2) * (-3. - 4. * C4 * s) +
              s * (3. + 2. * C4 * s)) -
        (0.25 * pow(-2 + delta, 2) * pow(mpion, 2) *
         (pow(mpion, 4) + pow(pow(mrho, 2) - t, 2) -
          2 * pow(mpion, 2) * (pow(mrho, 2) + t))) /
            (pow(mrho, 2) * pow(pow(mpion, 2) - t, 2)) +
        (pow(delta, 2) *
         (0.5 * pow(mpion, 4) * pow(mrho, 2) + 0.125 * pow(mrho, 6) +
          pow(mpion, 2) * (1.25 * pow(mrho, 4) - 0.375 * pow(s, 2) +
                           pow(mrho, 2) * (0.125 * s - 1. * t)) +
          pow(mrho, 4) * (-0.375 * s - 0.5 * t) +
          pow(s, 2) * (0.125 * s + 0.125 * t) +
          pow(mrho, 2) *
              (0.0625 * pow(s, 2) + 0.375 * s * t + 0.5 * pow(t, 2)))) /
            pow(mrho, 6) +
        (delta * (2. * C4 * pow(mrho, 6) +
                  pow(mpion, 2) * (3. * C4 * pow(mrho, 4) + 0.25 * s +
                                   pow(mrho, 2) * (-1.25 - 1. * C4 * s)) +
                  s * (-0.25 * s - 0.25 * t) +
                  pow(mrho, 4) * (-0.75 - 1.5 * C4 * s - 3. * C4 * t) +
                  pow(mrho, 2) *
                      (1.25 * t + s * (0.25 - 0.5 * C4 * s + 1. * C4 * t)))) /
            pow(mrho, 4) -
        (0.125 * (-2 + delta) * (eta1 - eta2) *
         (-(eta2 * (pow(mpion, 2) + t) *
            (pow(mpion, 4) + t * (-pow(mrho, 2) + 2 * s + t) -
             pow(mpion, 2) * (pow(mrho, 2) + 2 * t))) +
          eta1 * (2 * pow(mpion, 6) +
                  pow(mpion, 4) * (-2 * pow(mrho, 2) + s - 4 * t) +
                  s * t * (-pow(mrho, 2) + t) +
                  pow(mpion, 2) * (2 * pow(mrho, 4) + 2 * t * (s + t) -
                                   pow(mrho, 2) * (s + 2 * t))))) /
            ((-pow(ma1, 2) + t) * (-pow(mpion, 2) + t)) +
        (0.03125 * pow(eta1 - eta2, 2) *
         (-2 * eta1 * eta2 *
              (pow(mpion, 8) - 4 * pow(mpion, 6) * t +
               pow(t, 2) * (-pow(mrho, 4) - 2 * pow(mrho, 2) * s +
                            2 * pow(s, 2) + 2 * s * t + pow(t, 2)) -
               2 * pow(mpion, 2) * t *
                   (-2 * pow(mrho, 4) + pow(mrho, 2) * s + 2 * t * (s + t)) +
               pow(mpion, 4) * (-pow(mrho, 4) + 2 * t * (s + 3 * t))) +
          pow(eta2, 2) *
              (pow(mpion, 8) - 2 * pow(mpion, 6) * (pow(mrho, 2) + 2 * t) +
               pow(t, 2) * (pow(mrho, 4) + 2 * pow(s, 2) + 2 * s * t +
                            pow(t, 2) + 2 * pow(mrho, 2) * (-s + t)) -
               2 * pow(mpion, 2) * t *
                   (2 * t * (s + t) + pow(mrho, 2) * (s + 3 * t)) +
               pow(mpion, 4) * (pow(mrho, 4) + 6 * pow(mrho, 2) * t +
                                2 * t * (s + 3 * t))) +
          pow(eta1, 2) *
              (pow(mpion, 8) + 2 * pow(mpion, 6) * (pow(mrho, 2) - 2 * t) -
               2 * pow(mpion, 2) * (pow(mrho, 2) - s - t) *
                   (pow(mrho, 4) + pow(mrho, 2) * t - 2 * pow(t, 2)) +
               t * (-pow(mrho, 2) + t) *
                   (2 * pow(s, 2) + 2 * s * t + pow(t, 2) -
                    pow(mrho, 2) * (2 * s + t)) +
               pow(mpion, 4) * (pow(mrho, 4) - 2 * pow(mrho, 2) * (s + 3 * t) +
                                2 * t * (s + 3 * t))))) /
            pow(pow(ma1, 2) - t, 2) -
        (0.0625 * (eta1 - eta2) *
         (eta2 *
              (-4 * delta * pow(mpion, 6) +
               4 * pow(mpion, 4) *
                   (pow(mrho, 2) - 2 * C4 * pow(mrho, 4) + 3 * delta * t) +
               pow(mpion, 2) * (delta * (s - 6 * t) * (s + 2 * t) -
                                (2 + delta) * pow(mrho, 2) * (s + 4 * t) +
                                2 * pow(mrho, 4) * (-1 + delta + 8 * C4 * t)) +
               t * (-8 * C4 * pow(mrho, 6) +
                    pow(mrho, 4) * (6 - 4 * delta + 16 * C4 * s - 8 * C4 * t) +
                    pow(mrho, 2) * (-(s * (2 + delta + 8 * C4 * s)) +
                                    4 * (1 + delta) * t) +
                    delta * (3 * pow(s, 2) + 4 * s * t + 4 * pow(t, 2)))) +
          eta1 *
              (4 * delta * pow(mpion, 6) - 8 * C4 * pow(mrho, 6) * t +
               delta * (pow(s, 3) - 4 * pow(s, 2) * t - 6 * s * pow(t, 2) -
                        4 * pow(t, 3)) -
               2 * pow(mpion, 4) *
                   ((2 - 5 * delta) * pow(mrho, 2) - 4 * C4 * pow(mrho, 4) +
                    delta * (s + 6 * t)) +
               2 * pow(mrho, 4) *
                   (s - delta * s + t * (2 - delta + 4 * C4 * t)) +
               pow(mrho, 2) *
                   (8 * delta * s * t + 2 * (-2 + 3 * delta) * pow(t, 2) +
                    pow(s, 2) * (-2 + delta + 8 * C4 * t)) -
               2 * pow(mpion, 2) *
                   (-8 * C4 * pow(mrho, 6) + 2 * delta * (s - 3 * t) * (s + t) -
                    pow(mrho, 2) * ((2 + delta) * s + (4 - 8 * delta) * t) +
                    pow(mrho, 4) * (4 + 8 * C4 * (s + t)))))) /
            (pow(mrho, 2) * (-pow(ma1, 2) + t)) +
        (0.0625 * pow(-2 * pow(mrho, 2) + delta * s, 2) *
         (8 * pow(mpion, 4) * pow(mrho, 2) + 2 * pow(mrho, 6) + pow(s, 3) +
          8 * pow(mrho, 2) * t * (s + t) - pow(mrho, 4) * (7 * s + 8 * t) +
          4 * pow(mpion, 2) *
              (5 * pow(mrho, 4) - pow(s, 2) - 4 * pow(mrho, 2) * t)) *
         HeavisideTheta(-mrho + sqrt(s))) /
            (pow(mrho, 6) * pow(pow(mrho, 2) - s, 2)) -
        (0.0625 * (eta1 - eta2) * (2 * pow(mrho, 2) - delta * s) *
         (-(eta2 * (2 * pow(mpion, 2) + pow(mrho, 2) - s - 2 * t) *
            (2 * pow(mpion, 4) + pow(mpion, 2) * (-pow(mrho, 2) + s - 4 * t) +
             t * (3 * pow(mrho, 2) + s + 2 * t))) +
          eta1 * (4 * pow(mpion, 6) - pow(mrho, 4) * s + pow(s, 3) +
                  2 * pow(mpion, 4) * (5 * pow(mrho, 2) - s - 6 * t) -
                  2 * (pow(mrho, 4) - 4 * pow(mrho, 2) * s + pow(s, 2)) * t +
                  6 * (pow(mrho, 2) - s) * pow(t, 2) - 4 * pow(t, 3) -
                  4 * pow(mpion, 2) *
                      (4 * pow(mrho, 2) * t + (s - 3 * t) * (s + t)))) *
         HeavisideTheta(-mrho + sqrt(s))) /
            (pow(mrho, 2) * (pow(mrho, 2) - s) * (pow(ma1, 2) - t)) -
        (3. * (1. * pow(mrho, 2) - 0.5 * delta * s) *
         (delta *
              (0.666667 * pow(mpion, 4) * pow(mrho, 2) +
               0.166667 * pow(mrho, 6) +
               pow(mpion, 2) * (1.66667 * pow(mrho, 4) - 0.416667 * pow(s, 2) +
                                pow(mrho, 2) * (0.0833333 * s - 1.33333 * t)) +
               pow(mrho, 4) * (-0.541667 * s - 0.666667 * t) +
               pow(s, 2) * (0.125 * s + 0.0833333 * t) +
               pow(mrho, 2) * (-0.0833333 * pow(s, 2) + 0.583333 * s * t +
                               0.666667 * pow(t, 2))) +
          pow(mrho, 2) *
              (1. * C4 * pow(mrho, 6) +
               pow(mpion, 2) *
                   (2. * C4 * pow(mrho, 4) + 0.166667 * s +
                    pow(mrho, 2) * (-0.833333 - 0.666667 * C4 * s)) +
               s * (-0.0833333 * s - 0.166667 * t) +
               pow(mrho, 4) * (-0.416667 - 1.33333 * C4 * s - 2. * C4 * t) +
               pow(mrho, 2) * (0.833333 * t + s * (0.5 + 0.333333 * C4 * s +
                                                   0.666667 * C4 * t)))) *
         HeavisideTheta(-mrho + sqrt(s))) /
            (pow(mrho, 8) - 1. * pow(mrho, 6) * s) +
        ((0.125 * (-2 + delta) *
          (pow(mpion, 4) * ((-2 + 4 * delta) * pow(mrho, 2) +
                            8 * C4 * pow(mrho, 4) + 5 * delta * s) -
           8 * C4 * pow(mrho, 6) * t + delta * s * t * (s + t) +
           pow(mrho, 2) * (delta * s * (s - 3 * t) - 2 * t * (s + t)) +
           2 * pow(mrho, 4) *
               ((-1 + delta) * s + t + 4 * C4 * t * (2 * s + t)) -
           pow(mpion, 2) * (8 * C4 * pow(mrho, 6) + delta * s * (s + 6 * t) +
                            2 * pow(mrho, 4) * (-3 + 8 * C4 * t) +
                            pow(mrho, 2) * ((-2 + 9 * delta) * s +
                                            4 * (-1 + delta) * t)))) /
             (-pow(mpion, 2) + t) -
         (0.125 * (-2. + delta) * (-2. * pow(mrho, 2) + delta * s) *
          (pow(mpion, 4) * (4. * pow(mrho, 2) + 4. * s) +
           pow(mpion, 2) *
               (pow(mrho, 2) * (-7. * s - 4. * t) + s * (-1. * s - 4. * t)) +
           s * (pow(mrho, 4) + pow(mrho, 2) * (s - 1. * t) + s * t)) *
          HeavisideTheta(-mrho + sqrt(s))) /
             ((pow(mrho, 2) - 1. * s) * (pow(mpion, 2) - 1. * t))) /
            pow(mrho, 4))) /
      (16. * Pi * s * (-4 * pow(mpion, 2) + s));

  return diff_sigma / spin_deg_factor;
}
