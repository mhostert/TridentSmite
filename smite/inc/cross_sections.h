#ifndef CROSS_S_
#define CROSS_S_

#include <cmath>
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime> 
#include <cstdlib> // Need it for rand() and etc
#include <functional> // Need it for function as argument
#include <numeric>

#include "cuba.h"
#include "phase_space.h"


// functions for the integrands (using Czys et al formalism)
/* -------------------------------------------------------------------------------------- */

int CohFromIntToPhysical_Efixed(const cubareal xx[], void *MC);
int CohFromIntToPhysical(const cubareal xx[], void *MC);
int DifFromIntToPhysical_Efixed(const cubareal xx[], void *MC);
int DifFromIntToPhysical(const cubareal xx[], void *MC);

// Flux functions
/* -------------------------------------------------------------------------------------- */

long double hT_coh(int A, int Z, long double M, long double Enu, long double shat, long double Q2);
long double hL_coh(int A, int Z, long double M, long double Enu, long double shat, long double Q2);



// EPA virtual Cross sections for all channels
/* -------------------------------------------------------------------------------------- */

long double dsigma_dPS(int nu_alpha, int l1, int l2, int A, int Z, long double Enu, long double s, long double phi, long double theta, long double t, long double l, long double q, std::vector<long double> &BSM_params);
long double dsigma_dPS_diff(int nu_alpha, int l1, int l2, int A, int Z, long double Enu, long double s, long double phi, long double theta, long double t, long double l, long double q, std::vector<long double> &BSM_params);


// FORM FACTORS
/* -------------------------------------------------------------------------------------- */
long double FF_WS(long double q, long double A);
long double u1_WS(long double x1, long double A);
long double x1_WS(long double u1, long double A); 

long double H1_n(long double q);
long double H2_n(long double q);
long double H1_p(long double q);
long double H2_p(long double q);

long double F(long double q2, long double A);
long double F_diffractive(long double q2);

long double F_PAULI_BLOCK(long double x1);

#endif
