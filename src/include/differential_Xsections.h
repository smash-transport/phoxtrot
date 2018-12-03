// differential_Xsections.h

#ifndef DIFFERENTIAL_XSECTIONS_H_
#define DIFFERENTIAL_XSECTIONS_H_

#include <cmath>
#include <fstream>
#include <iostream>

//C11: (π + ρ0 -> (π, ρ, a1) -> π + γ)
double inline diff_xsection_C11(double t, double s);

//C12: (π0 + ρ -> (π, ρ, a1) -> π + γ)
double inline diff_xsection_C12(double t, double s);

//C13: (π + ρ -> (π, ρ, a1) -> π0 + γ)
double inline diff_xsection_C13(double t, double s);

//C14: (π0 + ρ0 -> (ω) -> π0 + γ)
double inline diff_xsection_C14(double t, double s);

//C15: (π + ρ -> (ω) -> π0 + γ)
double inline diff_xsection_C15(double t, double s);

//C16: (π0 + ρ -> (ω) -> π + γ)
double inline diff_xsection_C16(double t, double s);

//C21: (π + π -> (π, ρ, a1) -> ρ0 + γ)
double inline diff_xsection_C21(double t, double s);

//C22: (π + π0 -> (π, ρ, a1) -> ρ + γ)
double inline diff_xsection_C22(double t, double s);


#endif
