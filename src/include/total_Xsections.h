// total_Xsections.h

#ifndef TOTAL_XSECTIONS_H_
#define TOTAL_XSECTIONS_H_

#include <cmath>
#include <fstream>
#include <iostream>


//C11: (π + ρ0 -> (π, ρ, a1) -> π + γ)
double inline total_xsection_C11(double s);

//C12: (π0 + ρ -> (π, ρ, a1) -> π + γ)
double inline total_xsection_C11(double s);

//C13: (π + ρ -> (π, ρ, a1) -> π0 + γ)
double inline total_xsection_C11(double s);

//C14: (π0 + ρ0 -> (ω) -> π0 + γ)
double inline total_xsection_C11(double s);

//C15: (π + ρ -> (ω) -> π0 + γ)
double inline total_xsection_C11(double s);

//C16: (π0 + ρ -> (ω) -> π + γ)
double inline total_xsection_C11(double s);

//C21: (π + π -> (π, ρ, a1) -> ρ0 + γ)
double inline total_xsection_C11(double s);

//C22: (π + π0 -> (π, ρ, a1) -> ρ + γ)
double inline total_xsection_C11(double s);


#endif
