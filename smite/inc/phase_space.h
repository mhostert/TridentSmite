#ifndef PS_H_
#define PS_H_

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iomanip>      // std::setprecision

#include "constants.h"
#include "cuba.h"

std::vector<long double> phase_space_flux(const cubareal xx[], void *MC);
std::vector<long double> phase_space_Efixed(const cubareal xx[], void *MC);


#endif