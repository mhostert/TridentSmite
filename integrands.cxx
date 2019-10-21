#include "integrands.h"

int counter;

///////////////////////////////////////////////////////
int coh_BSM_E(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata,
  const int *number, const int core, const double *w, const int *iter) {
///////////////////////////////////////////////////////

  std::vector<long double> xphysical;

  // Get physical variables from integration variables
  xphysical = CohFromIntToPhysical_Efixed( xx, userdata);

  long double x1    = xphysical[0];
  long double x2    = xphysical[1];
  long double x3    = xphysical[2];
  long double x4    = xphysical[3];
  long double x5    = xphysical[4];
  long double x6    = xphysical[5];
  long double x7    = xphysical[6];
  long double x8    = xphysical[7];
  long double Enu    = xphysical[8];
  long double dsigma    = xphysical[9]; //GeV^-2

  // Changing units
  dsigma *= (GeV2_to_cm2*1e48);
 
  ///////////////////////////////////////////////////////////////////////////////////
  // CONVOLVING THE FLUX
  tridentMC* myMC = (tridentMC *)userdata;
  ///////////////////////////////////////////////////////////////////////////////////
  // final integrand
  ff[0] = dsigma;

  ///////////////////////////////////////////////////////////////////////////////////
  // print phase space point to file with integrand and weights
  myMC->my_samples.precision(17);
  // if (myMC->SAMPLES_FLAG == PRINT_SAMPLES 
        // && NNEW*(*iter)*(*iter +1)/2.0 > SKIP_PRINTS * MAXEVAL)
  if (myMC->SAMPLES_FLAG == PRINT_SAMPLES)  //       && counter > SKIP_PRINTS * MAXEVAL)
  {
   myMC->my_samples<<x1<<" "<<x2<<" "<<x3<<" "<<x4<<" "<<x5<<" "<<x6<<" "<<x7<<" "<<x8<<" "<<Enu<<" "<<ff[0]<<" "<<*w<<" "<<*iter<<std::endl;

  }

  // Two points per integrand call in Suave!! 
  // counter=1+counter;

  return 0;

}


///////////////////////////////////////////////////////
int coh_BSM_flux(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata,
  const int *number, const int core, const double *w, const int *iter) {
///////////////////////////////////////////////////////
  
  std::vector<long double> xphysical;

  // Get physical variables from integration variables
  xphysical = CohFromIntToPhysical( xx, userdata);

  long double x1    = xphysical[0];
  long double x2    = xphysical[1];
  long double x3    = xphysical[2];
  long double x4    = xphysical[3];
  long double x5    = xphysical[4];
  long double x6    = xphysical[5];
  long double x7    = xphysical[6];
  long double x8    = xphysical[7];
  long double Enu    = xphysical[8];
  long double dsigma    = xphysical[9]; //GeV^-2

  // Changing units
  dsigma *= (GeV2_to_cm2*1e48);
 
  ///////////////////////////////////////////////////////////////////////////////////
  // CONVOLVING THE FLUX
  tridentMC* myMC = (tridentMC *)userdata;
  dsigma *= (myMC->dPHIdE[floor(Enu/myMC->dE)-myMC->Ei]);
  dsigma *= (myMC->Emax - myMC->Emin); // jacobian of energy integration

  ///////////////////////////////////////////////////////////////////////////////////
  // final integrand
  ff[0] = dsigma;

  ///////////////////////////////////////////////////////////////////////////////////
  // print phase space point to file with integrand and weights
  myMC->my_samples.precision(17);
  // if (myMC->SAMPLES_FLAG == PRINT_SAMPLES 
        // && NNEW*(*iter)*(*iter +1)/2.0 > SKIP_PRINTS * MAXEVAL)
  if (myMC->SAMPLES_FLAG == PRINT_SAMPLES) //       && counter > SKIP_PRINTS * MAXEVAL)
  {
   myMC->my_samples<<x1<<" "<<x2<<" "<<x3<<" "<<x4<<" "<<x5<<" "<<x6<<" "<<x7<<" "<<x8<<" "<<Enu<<" "<<ff[0]<<" "<<*w<<" "<<*iter<<std::endl;
  }

  // Two points per integrand call in Suave!! 
  // counter=1+counter;
  // std::cout<<"counter: "<<counter<<std::endl;

  return 0;

}


///////////////////////////////////////////////////////

///////////////////////////////////////////////////////
int dif_BSM_E(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata,
  const int *number, const int core, const double *w, const int *iter) {
///////////////////////////////////////////////////////

  std::vector<long double> xphysical;

  // Get physical variables from integration variables
  xphysical = DifFromIntToPhysical_Efixed( xx, userdata);

  long double x1    = xphysical[0];
  long double x2    = xphysical[1];
  long double x3    = xphysical[2];
  long double x4    = xphysical[3];
  long double x5    = xphysical[4];
  long double x6    = xphysical[5];
  long double x7    = xphysical[6];
  long double x8    = xphysical[7];
  long double Enu    = xphysical[8];
  long double dsigma    = xphysical[9]; //GeV^-2

  // Changing units
  dsigma *= (GeV2_to_cm2*1e48);
 
  ///////////////////////////////////////////////////////////////////////////////////
  // CONVOLVING THE FLUX
  tridentMC* myMC = (tridentMC *)userdata;
  ///////////////////////////////////////////////////////////////////////////////////
  // final integrand
  ff[0] = dsigma;

  ///////////////////////////////////////////////////////////////////////////////////
  // print phase space point to file with integrand and weights
  myMC->my_samples.precision(17);
  // if (myMC->SAMPLES_FLAG == PRINT_SAMPLES 
        // && NNEW*(*iter)*(*iter +1)/2.0 > SKIP_PRINTS * MAXEVAL)
  if (myMC->SAMPLES_FLAG == PRINT_SAMPLES)  //       && counter > SKIP_PRINTS * MAXEVAL)
  {
    myMC->my_samples<<x1<<" "<<x2<<" "<<x3<<" "<<x4<<" "<<x5<<" "<<x6<<" "<<x7<<" "<<x8<<" "<<Enu<<" "<<ff[0]<<" "<<*w<<" "<<*iter<<std::endl;
  }

  // Two points per integrand call in Suave!! 
  // counter=1+counter;

  return 0;

}


///////////////////////////////////////////////////////
int dif_BSM_flux(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata,
  const int *number, const int core, const double *w, const int *iter) {
///////////////////////////////////////////////////////
  
  std::vector<long double> xphysical;

  // Get physical variables from integration variables
  xphysical = DifFromIntToPhysical( xx, userdata);

  long double x1    = xphysical[0];
  long double x2    = xphysical[1];
  long double x3    = xphysical[2];
  long double x4    = xphysical[3];
  long double x5    = xphysical[4];
  long double x6    = xphysical[5];
  long double x7    = xphysical[6];
  long double x8    = xphysical[7];
  long double Enu    = xphysical[8];
  long double dsigma = xphysical[9]; //GeV^-2

  // Changing units
  dsigma *= (GeV2_to_cm2*1e48);
 
  ///////////////////////////////////////////////////////////////////////////////////
  // CONVOLVING THE FLUX
  tridentMC* myMC = (tridentMC *)userdata;
  dsigma *= (myMC->dPHIdE[floor(Enu/myMC->dE)-myMC->Ei]);
  dsigma *= (myMC->Emax - myMC->Emin); // jacobian of energy integration

  ///////////////////////////////////////////////////////////////////////////////////
  // final integrand

  ff[0] = dsigma;

  ///////////////////////////////////////////////////////////////////////////////////
  // print phase space point to file with integrand and weights
  myMC->my_samples.precision(17);
  // if (myMC->SAMPLES_FLAG == PRINT_SAMPLES 
        // && NNEW*(*iter)*(*iter +1)/2.0 > SKIP_PRINTS * MAXEVAL)
  if (myMC->SAMPLES_FLAG == PRINT_SAMPLES)  //       && counter > SKIP_PRINTS * MAXEVAL)
  {
   myMC->my_samples<<x1<<" "<<x2<<" "<<x3<<" "<<x4<<" "<<x5<<" "<<x6<<" "<<x7<<" "<<x8<<" "<<Enu<<" "<<ff[0]<<" "<<*w<<" "<<*iter<<std::endl;

  }

  // Two points per integrand call in Suave!! 
  // counter=1+counter;

  return 0;

}





///////////////////////////////////////////////////////
int Coherent_EPA_E(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata,
  const int *number, const int core, const double *w, const int *iter) {
///////////////////////////////////////////////////////

  long double s     = xx[0];
  long double phi   = xx[1];
  long double theta = xx[2];
  long double t     = xx[3];
  long double l     = xx[4];
  long double q     = xx[5];

  tridentMC* myMC = (tridentMC *)userdata;
  
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

  // std::cout << "dsigma = "<< ff[0] << ",  jacob = " << Jacob <<std::endl;
  ff[0] *= Jacob * 1e48;


 // std::cout<< <<std::endl;

  // myMC->my_samples         <<s_true<<" "
  //                            <<phi_true<<" "
  //                            <<theta_true<<" "
  //                            <<t_true<<" "
  //                            <<l_true+t_true/*-SQR(Mzprime)*exp(l_true)+t_true+SQR(exp(q_true))+SQR(Mzprime)*/<<" " //////////////////// VARIABLE SUBSTITUTION !!!!!
  //                            <<q_true/*exp(q_true)*/<<" "
  //                            <<Enu_true<<" "
  //                            <<ff[0]<<" "
  //                            <<*w<<" "
  //                            <<*iter
  //                            <<std::endl;


  return 0;
}
