#include "integrands.h"

///////////////////////////////////////////////////////
int coh_BSM_E(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *MC,
  const int *number, const int core, const double *w, const int *iter) {
///////////////////////////////////////////////////////

  tridentMC* myMC = (tridentMC *)MC;

  // Get physical variables from integration variables
  CohFromIntToPhysical_Efixed( xx, MC);

  long double x1    = myMC->xphys[0];
  long double x2    = myMC->xphys[1];
  long double x3    = myMC->xphys[2];
  long double x4    = myMC->xphys[3];
  long double x5    = myMC->xphys[4];
  long double x6    = myMC->xphys[5];
  long double x7    = myMC->xphys[6];
  long double x8    = myMC->xphys[7];
  long double Enu    = myMC->xphys[8];
  long double dsigma    = myMC->xphys[9]; // units of yb = yoctobarn = 1e-48 cm^2


  ///////////////////////////////////////////////////////////////////////////////////
  // final integrand
  ff[0] = dsigma;


  // save samples and some info about event weight and iteration number
  myMC->xphys[10]=(*w);
  myMC->xphys[11]=(*iter);
  myMC->cuba_samples.push_back(myMC->xphys);

  return 0;

}


///////////////////////////////////////////////////////
int coh_BSM_flux(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *MC,
  const int *number, const int core, const double *w, const int *iter) {
///////////////////////////////////////////////////////

  tridentMC* myMC = (tridentMC *)MC;

  // Get physical variables from integration variables
  CohFromIntToPhysical( xx, MC);

  long double x1    = myMC->xphys[0];
  long double x2    = myMC->xphys[1];
  long double x3    = myMC->xphys[2];
  long double x4    = myMC->xphys[3];
  long double x5    = myMC->xphys[4];
  long double x6    = myMC->xphys[5];
  long double x7    = myMC->xphys[6];
  long double x8    = myMC->xphys[7];
  long double Enu    = myMC->xphys[8];
  long double dsigma    = myMC->xphys[9]; // units of yb = yoctobarn = 1e-48 cm^2
 
  ///////////////////////////////////////////////////////////////////////////////////
  // CONVOLVING THE FLUX
  dsigma *= (myMC->dPHIdE[floor(Enu/myMC->dE)-myMC->Ei]);
  dsigma *= (myMC->Emax - myMC->Emin); // jacobian of energy integration

  ///////////////////////////////////////////////////////////////////////////////////
  // final integrand
  ff[0] = dsigma;

  myMC->xphys[10]=(*w);
  myMC->xphys[11]=(*iter);
  myMC->cuba_samples.push_back(myMC->xphys);

  return 0;

}


///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
int dif_BSM_E(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *MC,
  const int *number, const int core, const double *w, const int *iter) {
///////////////////////////////////////////////////////

  tridentMC* myMC = (tridentMC *)MC;

  // Get physical variables from integration variables
  DifFromIntToPhysical_Efixed( xx, MC);


  long double x1    = myMC->xphys[0];
  long double x2    = myMC->xphys[1];
  long double x3    = myMC->xphys[2];
  long double x4    = myMC->xphys[3];
  long double x5    = myMC->xphys[4];
  long double x6    = myMC->xphys[5];
  long double x7    = myMC->xphys[6];
  long double x8    = myMC->xphys[7];
  long double Enu    = myMC->xphys[8];
  long double dsigma    = myMC->xphys[9]; // units of yb = yoctobarn = 1e-48 cm^2


  ///////////////////////////////////////////////////////////////////////////////////
  // final integrand
  ff[0] = dsigma;


  // save samples and some info about event weight and iteration number
  myMC->xphys[10]=(*w);
  myMC->xphys[11]=(*iter);
  myMC->cuba_samples.push_back(myMC->xphys);

  return 0;

}


///////////////////////////////////////////////////////
int dif_BSM_flux(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *MC,
  const int *number, const int core, const double *w, const int *iter) {
///////////////////////////////////////////////////////

  tridentMC* myMC = (tridentMC *)MC;

  // Get physical variables from integration variables
  DifFromIntToPhysical( xx, MC);


  long double x1    = myMC->xphys[0];
  long double x2    = myMC->xphys[1];
  long double x3    = myMC->xphys[2];
  long double x4    = myMC->xphys[3];
  long double x5    = myMC->xphys[4];
  long double x6    = myMC->xphys[5];
  long double x7    = myMC->xphys[6];
  long double x8    = myMC->xphys[7];
  long double Enu    = myMC->xphys[8];
  long double dsigma    = myMC->xphys[9]; // units of yb = yoctobarn = 1e-48 cm^2


  ///////////////////////////////////////////////////////////////////////////////////
  // CONVOLVING THE FLUX
  dsigma *= (myMC->dPHIdE[floor(Enu/myMC->dE)-myMC->Ei]);
  dsigma *= (myMC->Emax - myMC->Emin); // jacobian of energy integration

  ///////////////////////////////////////////////////////////////////////////////////
  // final integrand
  ff[0] = dsigma;

  // save samples and some info about event weight and iteration number
  myMC->xphys[10]=(*w);
  myMC->xphys[11]=(*iter);
  myMC->cuba_samples.push_back(myMC->xphys);

  return 0;

}





///////////////////////////////////////////////////////
int Coherent_EPA_E(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *MC,
  const int *number, const int core, const double *w, const int *iter) {
///////////////////////////////////////////////////////

  long double s     = xx[0];
  long double phi   = xx[1];
  long double theta = xx[2];
  long double t     = xx[3];
  long double l     = xx[4];
  long double q     = xx[5];

  tridentMC* myMC = (tridentMC *)MC;
  
  long double ml1 = myMC->ml1;
  long double ml2 = myMC->ml2;
  long double A = myMC->A;
  long double Z = myMC->Z;
  
  long double Enu_true = myMC->nu_energy;

  // This defines the cuts on Q2 by using the new physics parameter vector
  long double Qmax_coh;
  Qmax_coh = myMC->Qcut_coh;
  
  ///////////////////////////////////////////
  // Variable substitution to unit hypercube  
  long double q_true_min = SQR(ml1+ml2)/2.0/Enu_true;

  long double q_true = q*(Qmax_coh - q_true_min) + q_true_min;
  long double s_true = s*(2*Enu_true*q_true - SQR(ml1+ml2))+ SQR(ml1+ml2);
  long double phi_true = phi*2*M_PI;
  long double theta_true = theta*M_PI;
  long double t_true = t*(s_true - SQR(ml1+ml2)) + SQR(ml1+ml2); //VARIABLE SUBSTITUTION
  long double l_true = l*(t_true - SQR(ml1+ml2)) + SQR(ml1+ml2)-t_true;


  long double Jacob;

 Jacob = (M_PI) 
         * (2*M_PI) 
           * (2*Enu_true*q_true - SQR(ml1+ml2)) 
             * (t_true - SQR(ml1+ml2))  //VARIABLE SUBSTITUTION
               * (s_true-SQR(ml1+ml2))  //VARIABLE SUBSTITUTION
                 * (Qmax_coh - q_true_min);

 
  ff[0] = dsigma_dPS( myMC->nu_alpha, 
                      myMC->l1, 
                      myMC->l2, 
                      A,
                      Z,
                      Enu_true,
                      s_true,
                      phi_true,
                      theta_true,
                      t_true,
                      l_true+t_true, //VARIABLE SUBSTITUTION
                      // -SQR(Mzprime)*exp(l_true)+t_true+SQR(exp(q_true))+SQR(Mzprime),      
                      q_true,
                      //exp(q_true),
                      myMC->params);

  ff[0] *= Jacob * 1e48;

  return 0;
}
