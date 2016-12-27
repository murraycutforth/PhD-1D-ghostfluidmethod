#ifndef MISC_H
#define MISC_H


#include "eos.hpp"
#include <blitz/array.h>


blitz::Array<double,1> euler_flux (double rho, double u, double P, double E);

blitz::Array<double,1> euler_flux (blitz::Array<double,1> cv, std::shared_ptr<eos_base> eos);

double specific_ie_cv (blitz::Array<double,1> state);

bool is_state_physical (blitz::Array<double,1> state);

double gaussian_function (double A, double mu, double sigma, double x);


#endif
