#ifndef INTEGRANDS_H_
#define INTEGRANDS_H_

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>
#include <sstream>
#include <ctime>

#include "cuba.h"
#include "constants.h"
#include "cross_sections.h"
#include "observables.h"

///////////////////////////////////////////////////////
int coh_BSM_E(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata,
  const int *number, const int core, const double *w, const int *iter);
///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
int coh_BSM_flux(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata,
  const int *number, const int core, const double *w, const int *iter);
///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
int dif_BSM_E(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata,
  const int *number, const int core, const double *w, const int *iter);
///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
int dif_BSM_flux(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata,
  const int *number, const int core, const double *w, const int *iter);
///////////////////////////////////////////////////////


///////////////////////////////////////////////////////
int Coherent_EPA_E(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata,
  const int *number, const int core, const double *w, const int *iter);
///////////////////////////////////////////////////////

#endif