// total_Xsections.h

#ifndef TWOTOTWO_TOTAL_XSECTIONS_H_
#define TWOTOTWO_TOTAL_XSECTIONS_H_

#include <cmath>
#include <fstream>
#include <iostream>

class TwoToTwo_Tot_Xsections {
  public:
    // C11: (π + ρ0 -> (π, ρ, a1) -> π + γ)
    static double total_xsection_C11(double s);

    // C12: (π0 + ρ -> (π, ρ, a1) -> π + γ)
    static double total_xsection_C12(double s);

    // C13: (π + ρ -> (π, ρ, a1) -> π0 + γ)
    static double total_xsection_C13(double s);

    // C14: (π0 + ρ0 -> (ω) -> π0 + γ)
    static double total_xsection_C14(double s);

    // C15: (π + ρ -> (ω) -> π0 + γ)
    static double total_xsection_C15(double s);

    // C16: (π0 + ρ -> (ω) -> π + γ)
    static double total_xsection_C16(double s);

    // C21: (π + π -> (π, ρ, a1) -> ρ0 + γ)
    static double total_xsection_C21(double s);

    // C22: (π + π0 -> (π, ρ, a1) -> ρ + γ)
    static double total_xsection_C22(double s);
};

#endif
