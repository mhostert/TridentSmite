#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

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
#include <iomanip>      // std::setprecision

#include "integrands.h"


//////////////////////////////////////////////////////
// Class that will contain all the userdata to be passed to INTEGRAND
class tridentMC{
  public:

    int nu_alpha;
    int l1,l2,A,Z;

    int IS_NUBAR = 0;

    long double ml1,ml2,Mn;

    // Will contain BSM parameters
    std::vector<long double> params;

    // Will contain the choice for what terms to include
    std::vector<long double> terms;

    // neutrino energy when not integrating over the flux
    long double nu_energy;

    // boundaries of Q integration
    long double Qcut_coh, Qmin_dif, Qmax_dif;

    // Channels axial and vector couplings in the SM
    long double Aijk, Vijk;
    long double CHARGE, gprimeA, gprimeV;

    std::vector<long double> dsigma;

    std::ofstream my_samples;
    std::ofstream my_sums;
    std::ifstream my_flux;

    int E_FLAG;
    int SAMPLES_FLAG;
    int TOTAL_DIAGRAMS;
    int PAULI_BLOCKING = W_BLOCKING;

    // Flux variables
    std::vector<long double> Evec, dPHIdE;
    long double dE, Emax, Emin;
    int Ei;
    
    // Vegas helper variables
    int comp, regions, neval, fail, counter;

    cubareal integral[NCOMP], error[NCOMP], chi2prob[NCOMP];
    cubareal integral_L[NCOMP], error_L[NCOMP], chi2prob_L[NCOMP];
    cubareal integral_T[NCOMP], error_T[NCOMP], chi2prob_T[NCOMP];

    tridentMC(int , int , int , long double );

    void open_flux_file(std::string , long double, long double);

    long double integrate_wflux(void *, int );
    long double integrate_wflux_wsamples(void *, int , std::string );

    long double integrate_energy(void *, long double  ,int );
    long double integrate_energy_wsamples(void *, long double  ,int , std::string );

    long double compute_total_xsec_energy(void * integrando, long double Enu, int ndim);

    long double final_xsec(long double gprimeV, long double gprimeA, long double CHARGE);    

    void generate_events(std::string , std::string );

    std::vector<std::vector<long double>> histogram(long double xmin,long double xmax, int NBINS, int VAR_FLAG, std::string eventsfile);


    void remove_events(std::string eventsfile);

};

#endif