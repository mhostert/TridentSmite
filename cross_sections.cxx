#include "cross_sections.h"
#include "integrator.h"


int CohFromIntToPhysical_Efixed(const cubareal xx[], void *MC){

  long double u1_s    = xx[0];
  long double u2_s    = xx[1];
  long double u3_s    = xx[2];
  long double PHI2_s  = xx[3];
  long double x5_s    = xx[4];
  long double u6_s    = xx[5];
  long double x7_s    = xx[6];
  long double x8_s    = xx[7];

  // Random angles needed to specify 4 vectors
  long double x7 = x7_s*2*M_PI;
  long double x8 = x8_s*2*M_PI;

  //////////////////////////////////////////////////
  // get data from user

  tridentMC* myMC = (tridentMC *)MC;
  
  long double ml1 = myMC->ml1;
  long double ml2 = myMC->ml2;
  long double A = myMC->A;
  long double Z = myMC->Z;
  long double Mn = myMC->Mn;
  long double mzprime = myMC->mzprime;


  long double Diag11 = myMC->terms[0];
  long double Diag22 = myMC->terms[1];
  long double Diag12 = myMC->terms[2];
  long double V2     = myMC->terms[3];
  long double A2     = myMC->terms[4];
  long double VA     = myMC->terms[5];
  long double BSM    = myMC->terms[6];

  ////////////////////////////////////////////////////
  // FIXED neutrino energy 
  long double Enu = myMC->nu_energy;

  //////////////////////////////////
  long double x1_u = 2*SQR(Enu) / ( 1 + 2 * Enu/Mn ) * ( 1 - SQR(ml1+ml2)/2.0/SQR(Enu) * (1 + Enu/Mn) +
                  sqrt( SQR(1 - SQR(ml1 + ml2)/2.0/SQR(Enu)*(1+Enu/Mn) ) 
                    - pow(ml1 + ml2, 4)/4.0/pow(Enu,4)*(1 + 2*Enu/Mn) )  );

  long double x1_l = pow(ml1+ml2,4)/(x1_u)/(1.0+2*Enu/Mn);

  long double u1_l = log(x1_l); 
  long double u1_u = log(x1_u); 
  
  long double u1 = u1_s*(u1_u-u1_l) + u1_l;

  long double x1 = exp(u1);
  //////////////////////////////////
  long double x2_l = 1.0/2.0 * ( x1 + SQR(ml1 + ml2) );
  long double x2_u = Enu * (sqrt( x1 + SQR(x1)/4.0/SQR(Mn) ) - x1/2.0/Mn );
  // long double x2 = x2_s*(x2_u-x2_l) + x2_l;

  long double u2_u = (1.0+2.0*Enu/Mn)*(x1_u-x1)*(x1-x1_l)  / ( SQR(ml1 + ml2) + x1*(1+Enu/Mn) + 2*Enu*sqrt(x1+SQR(x1)/4.0/SQR(Mn)) );
  long double u2_l = 0.0;

  long double u2 = u2_s*(u2_u-u2_l) + u2_l;

  long double x2 = 0.5*(u2 + SQR(ml1 + ml2) + x1);

  /////////////////
  long double x5_l = x2/2.0/(2*x2-x1) * (  2*x2 - x1 + SQR(ml2) - SQR(ml1)  
        - sqrt( SQR( 2*x2-x1 + SQR(ml2) - SQR(ml1) ) - 4*SQR(ml2)*(2*x2-x1) ) );
  
  long double x5_u = x2/2.0/(2*x2-x1) * (  2*x2 - x1 + SQR(ml2) - SQR(ml1)  
        + sqrt( SQR( 2*x2-x1 + SQR(ml2) - SQR(ml1) ) - 4*SQR(ml2)*(2*x2-x1) ) );
  

  long double x5 = x5_s*(x5_u-x5_l) + x5_l; 

  ////////////////
  long double x3_l = SQR(ml2)/2.0 * x2/x5 - x1/2.0 * x5/x2;

  long double x3_u = 0.5*(2*x2 - x1 + SQR(ml2) - SQR(ml1) )  - x5;
  
  long double u3_l = log(2*x3_l + x1);
  long double u3_u = log(2*x3_u + x1);

  long double u3 = u3_s*(u3_u-u3_l) + u3_l;

  long double x3 = 0.5*(exp(u3)- x1);

  ////////////////


  ////////////////
  // Useful Definitions in Frame p1vec + qvec - p3vec = 0
  ////////////////
  long double Wc2 = 2*(x2-x3-x5) - x1 + SQR(ml2);

  long double E1 =  (x2 - x5)/sqrt(Wc2) ;
  long double E4 =  (Wc2+SQR(ml1))/2.0/sqrt(Wc2) ;
  long double q0 =  (x2-x1-x3)/sqrt(Wc2) ;
  long double qvec = sqrt(SQR(q0) + x1)  ;
  long double p4vec = (Wc2 - SQR(ml1))/2.0/sqrt(Wc2)  ;

  long double Cq = (q0*E1 - x2)/qvec/E1;
  long double Sq = sqrt(1 - SQR(Cq));

  long double E2 = (Wc2 - SQR(ml1))/2/sqrt(Wc2);


  ////////////////

  long double m6, u6, u6_l, u6_u, prop, jacob_u6;
  long double m6_l = 0;
  long double m6_u = 2*E1*E2;

  switch ((int) BSM) 
  {
    case (SMonly):
      u6_l = m6_l;
      u6_u = m6_u;
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6 = u6;
      jacob_u6 = 1.0;

      V2     = myMC->terms[3];
      A2     = myMC->terms[4];
      VA     = myMC->terms[5];
      break;

    case (INTERFERENCE):

      //////////////////////////////////////////////////////////
      // BSM PROPAGATOR 
      // long double prop = (1.0/(-2*(m6/mzprime/mzprime) - 1.0))/mzprime/mzprime;
      u6_l = log(1.0/(2.0*m6_u + mzprime*mzprime));
      u6_u = log(1.0/(2.0*m6_l + mzprime*mzprime));
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6   = -mzprime*mzprime/2.0 + exp(-u6)/2.0;

      prop = -exp(u6);
      jacob_u6 = exp(-u6)/2.0;

      V2     = myMC->terms[3] * prop;
      A2     = myMC->terms[4] * prop;
      VA     = myMC->terms[5] * prop;
      break;

    case (BSMonly):

      u6_l = 1.0/(2.0*m6_u + mzprime*mzprime);
      u6_u = 1.0/(2.0*m6_l + mzprime*mzprime);
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6   = -mzprime*mzprime/2.0 + 1.0/2.0/u6;

      prop = -u6;
      jacob_u6 = 1.0/(u6*u6)/2.0;

      V2     = myMC->terms[3] * prop*prop;
      A2     = myMC->terms[4] * prop*prop;
      VA     = myMC->terms[5] * prop*prop;
      break;
    
    case (SMandBSM):
     
      u6_l = m6_l;
      u6_u = m6_u;
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6 = u6;
      jacob_u6 = 1.0;

      prop = 1.0/(2.0*m6 + mzprime*mzprime);

      V2     = (myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop) ;
      A2     = (myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop) ;
      VA     = (myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop);
      break;

    default:
      std::cout << "Could not determine what contributions (BSM) to include." << std::endl;
      break;
  }

  long double C2 = 1.0 - m6 / E1 / E2;
  long double S2 = sqrt(1.0 - SQR(C2));


  ////////////////
  long double PHI2 = PHI2_s*(2*M_PI-0.0) + 0.0;

  long double x4 = q0*E4 + (Sq*S2*cos(PHI2) + Cq*C2)*E2*qvec;

  long double C4 = (E4*q0 - x4)/(p4vec*qvec);
  long double S4 = sqrt(1.0 - SQR(C4));

  long double x6 = E1*E4 + E1*E2*C2; 

  //////////////////////////////////////////////////////////
  // Multiply by the appropriate Jacobians for our invariants and phase space factors
  long double Jacob =  (Wc2 - SQR(ml1))/8.0/Wc2/E1/E2/x2/pow(2*M_PI,6.0);

  // Jacobian due to change from CSW(+Matheus) to better integration variables AND VEGAS
  Jacob *=  (u1_u-u1_l)*exp(u1)*
            (u2_u-u2_l)*0.5*
            (u3_u-u3_l)*0.5*exp(u3)*
            (2*M_PI)*
            (x5_u-x5_l)*
            (u6_u-u6_l)*jacob_u6*
            (2*M_PI)*
            (2*M_PI);


  //////////////////////////////////////////////////////////
  // MATRIX ELEMENT ITSELF
  long double FormFactor = FF_WS(sqrt(x1), A);

  // already includes Z*Z for the coherent enhancement
  long double dsigma = (8*(alphaQED*alphaQED)*(FormFactor*FormFactor)*(Gf*Gf)*(Diag22*((x1 + 2*x4)*(x1 + 2*x4))*(A2*(-2*(ml2*ml2)*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) - 2*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) + 2*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5)*x5 - x1*x1*(x2*x2)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 3*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 4*ml1*(ml2*ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(x2 - x5 - x6) + 2*ml1*ml2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6) - 4*ml1*ml2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(x2 - x5 - x6) - 2*ml1*ml2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*(x2 - x5 - x6) - 4*ml1*ml2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5)*(x2 - x5 - x6) - 4*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x6 - 4*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(-x1 + x2 - x3 - x4)*x6 - 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*x6 + 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 + x1*x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6)*x6 + 2*Mn*(x1*x1)*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x3*x5*(x2 - x5 - x6)*x6 + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x3*x5*(x2 - x5 - x6)*x6 - 2*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + x1*x1*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) + V2*(-2*(ml2*ml2)*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) - 2*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) + 2*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5)*x5 - x1*x1*(x2*x2)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 3*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 4*ml1*(ml2*ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(x2 - x5 - x6) - 2*ml1*ml2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6) + 4*ml1*ml2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(x2 - x5 - x6) + 2*ml1*ml2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*(x2 - x5 - x6) + 4*ml1*ml2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5)*(x2 - x5 - x6) - 4*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x6 - 4*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(-x1 + x2 - x3 - x4)*x6 - 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*x6 + 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 + x1*x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6)*x6 + 2*Mn*(x1*x1)*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x3*x5*(x2 - x5 - x6)*x6 + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x3*x5*(x2 - x5 - x6)*x6 - 2*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + x1*x1*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) + 2*VA*(-2*(ml2*ml2)*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) - 2*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) + 2*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5)*x5 - x1*x1*(x2*x2)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 3*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 4*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x6 + 4*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(-x1 + x2 - x3 - x4)*x6 + 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*x6 - 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 - 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 - 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 - x1*x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6)*x6 - 2*Mn*(x1*x1)*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(x2 - x5 - x6)*x6 - 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x3*x5*(x2 - x5 - x6)*x6 - 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x3*x5*(x2 - x5 - x6)*x6 + Mn*x1*(x2*x2)*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - x1*x1*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6))) + Diag11*((x1 + 2*x3)*(x1 + 2*x3))*(A2*(-4*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x5 - 4*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x4*x5 - 2*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 2*(ml1*ml1)*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) - 2*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x4*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) - 4*(ml1*ml1*ml1)*ml2*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x4*(x2 - x5 - x6) - 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*x6 + 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 + x1*x1*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 - 2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 + 2*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6)*x6 + 2*ml1*ml2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x2 - x5 - x6)*x6 - 4*ml1*ml2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6)*x6 + x1*x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6)*x6 + 2*Mn*(x1*x1)*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x4*x5*(x2 - x5 - x6)*x6 + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x4*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) - 2*ml1*ml2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x2 - x5 - x6)*(x6*x6) - 4*ml1*ml2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x2 - x5 - x6)*(x6*x6) - x1*x1*(x2*x2)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 3*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) + V2*(-4*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x5 - 4*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x4*x5 - 2*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 2*(ml1*ml1)*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) - 2*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x4*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) + 4*(ml1*ml1*ml1)*ml2*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x4*(x2 - x5 - x6) - 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*x6 + 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 + x1*x1*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 - 2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 + 2*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6)*x6 - 2*ml1*ml2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x2 - x5 - x6)*x6 + 4*ml1*ml2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6)*x6 + x1*x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6)*x6 + 2*Mn*(x1*x1)*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x4*x5*(x2 - x5 - x6)*x6 + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x4*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) + 2*ml1*ml2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x2 - x5 - x6)*(x6*x6) + 4*ml1*ml2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x2 - x5 - x6)*(x6*x6) - x1*x1*(x2*x2)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 3*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) + 2*VA*(-4*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x5 - 4*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x4*x5 - 2*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*(ml1*ml1)*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) + 2*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x4*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) - 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*x6 + 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 + x1*x1*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 - 2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 - 2*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6)*x6 + x1*x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6)*x6 + 2*Mn*(x1*x1)*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x4*x5*(x2 - x5 - x6)*x6 + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x4*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) + x1*x1*(x2*x2)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 3*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6))) + Diag12*(x1 + 2*x3)*(x1 + 2*x4)*(2*(Enu*Enu)*(Mn*Mn)*x1*(-2*V2*x1*(x2*x2*x2*x2) + 4*V2*(x2*x2*x2*x2*x2) + V2*x1*(x2*x2*x2)*x3 + 2*VA*x1*(x2*x2*x2)*x3 - 6*V2*(x2*x2*x2*x2)*x3 - 4*VA*(x2*x2*x2*x2)*x3 + 2*V2*(x2*x2*x2)*(x3*x3) + 4*VA*(x2*x2*x2)*(x3*x3) + V2*x1*(x2*x2*x2)*x4 - 2*VA*x1*(x2*x2*x2)*x4 - 6*V2*(x2*x2*x2*x2)*x4 + 4*VA*(x2*x2*x2*x2)*x4 + 4*V2*(x2*x2*x2)*x3*x4 + 2*V2*(x2*x2*x2)*(x4*x4) - 4*VA*(x2*x2*x2)*(x4*x4) + 4*V2*(x1*x1)*(x2*x2)*x5 + 8*VA*(x1*x1)*(x2*x2)*x5 - 6*V2*x1*(x2*x2*x2)*x5 - 12*VA*x1*(x2*x2*x2)*x5 - 8*V2*(x2*x2*x2*x2)*x5 - 8*VA*(x2*x2*x2*x2)*x5 + 8*V2*x1*(x2*x2)*x3*x5 + 16*VA*x1*(x2*x2)*x3*x5 + 6*V2*(x2*x2*x2)*x3*x5 + 12*VA*(x2*x2*x2)*x3*x5 + 8*V2*x1*(x2*x2)*x4*x5 + 16*VA*x1*(x2*x2)*x4*x5 + 6*V2*(x2*x2*x2)*x4*x5 - 4*VA*(x2*x2*x2)*x4*x5 + 2*V2*(x2*x2)*x3*x4*x5 + 4*VA*(x2*x2)*x3*x4*x5 - 2*V2*(x2*x2)*(x4*x4)*x5 + 4*VA*(x2*x2)*(x4*x4)*x5 - 3*V2*(x1*x1)*x2*(x5*x5) - 6*VA*(x1*x1)*x2*(x5*x5) + 12*V2*x1*(x2*x2)*(x5*x5) + 24*VA*x1*(x2*x2)*(x5*x5) + 8*V2*(x2*x2*x2)*(x5*x5) + 16*VA*(x2*x2*x2)*(x5*x5) - 6*V2*x1*x2*x3*(x5*x5) - 12*VA*x1*x2*x3*(x5*x5) - 4*V2*(x2*x2)*x3*(x5*x5) - 8*VA*(x2*x2)*x3*(x5*x5) - 12*V2*x1*x2*x4*(x5*x5) - 24*VA*x1*x2*x4*(x5*x5) - 6*V2*x1*x2*(x5*x5*x5) - 12*VA*x1*x2*(x5*x5*x5) - 4*V2*(x2*x2)*(x5*x5*x5) - 8*VA*(x2*x2)*(x5*x5*x5) + 6*V2*x1*x4*(x5*x5*x5) + 12*VA*x1*x4*(x5*x5*x5) + 4*V2*(x1*x1)*(x2*x2)*x6 - 8*VA*(x1*x1)*(x2*x2)*x6 - 6*V2*x1*(x2*x2*x2)*x6 + 12*VA*x1*(x2*x2*x2)*x6 - 8*V2*(x2*x2*x2*x2)*x6 + 8*VA*(x2*x2*x2*x2)*x6 + 8*V2*x1*(x2*x2)*x3*x6 - 16*VA*x1*(x2*x2)*x3*x6 + 6*V2*(x2*x2*x2)*x3*x6 + 4*VA*(x2*x2*x2)*x3*x6 - 2*V2*(x2*x2)*(x3*x3)*x6 - 4*VA*(x2*x2)*(x3*x3)*x6 + 8*V2*x1*(x2*x2)*x4*x6 - 16*VA*x1*(x2*x2)*x4*x6 + 6*V2*(x2*x2*x2)*x4*x6 - 12*VA*(x2*x2*x2)*x4*x6 + 2*V2*(x2*x2)*x3*x4*x6 - 4*VA*(x2*x2)*x3*x4*x6 + 16*V2*x1*(x2*x2)*x5*x6 + 8*V2*(x2*x2*x2)*x5*x6 - 6*V2*x1*x2*x3*x5*x6 - 12*VA*x1*x2*x3*x5*x6 - 4*V2*(x2*x2)*x3*x5*x6 - 8*VA*(x2*x2)*x3*x5*x6 - 6*V2*x1*x2*x4*x5*x6 + 12*VA*x1*x2*x4*x5*x6 - 4*V2*(x2*x2)*x4*x5*x6 + 8*VA*(x2*x2)*x4*x5*x6 - 18*V2*x1*x2*(x5*x5)*x6 - 36*VA*x1*x2*(x5*x5)*x6 - 4*V2*(x2*x2)*(x5*x5)*x6 - 8*VA*(x2*x2)*(x5*x5)*x6 + 6*V2*x1*x3*(x5*x5)*x6 + 12*VA*x1*x3*(x5*x5)*x6 + 12*V2*x1*(x5*x5*x5)*x6 + 24*VA*x1*(x5*x5*x5)*x6 - 3*V2*(x1*x1)*x2*(x6*x6) + 6*VA*(x1*x1)*x2*(x6*x6) + 12*V2*x1*(x2*x2)*(x6*x6) - 24*VA*x1*(x2*x2)*(x6*x6) + 8*V2*(x2*x2*x2)*(x6*x6) - 16*VA*(x2*x2*x2)*(x6*x6) - 12*V2*x1*x2*x3*(x6*x6) + 24*VA*x1*x2*x3*(x6*x6) - 6*V2*x1*x2*x4*(x6*x6) + 12*VA*x1*x2*x4*(x6*x6) - 4*V2*(x2*x2)*x4*(x6*x6) + 8*VA*(x2*x2)*x4*(x6*x6) - 18*V2*x1*x2*x5*(x6*x6) + 36*VA*x1*x2*x5*(x6*x6) - 4*V2*(x2*x2)*x5*(x6*x6) + 8*VA*(x2*x2)*x5*(x6*x6) + 6*V2*x1*x4*x5*(x6*x6) - 12*VA*x1*x4*x5*(x6*x6) - 6*V2*x1*x2*(x6*x6*x6) + 12*VA*x1*x2*(x6*x6*x6) - 4*V2*(x2*x2)*(x6*x6*x6) + 8*VA*(x2*x2)*(x6*x6*x6) + 6*V2*x1*x3*(x6*x6*x6) - 12*VA*x1*x3*(x6*x6*x6) + 12*V2*x1*x5*(x6*x6*x6) - 24*VA*x1*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1)*ml2*V2*(x2*x2)*(-x2 + x5 + x6) + ml2*ml2*ml2*ml2*(x2*x2)*(2*VA*(x2 - x5 - x6) + V2*(-x5 + x6)) + ml1*ml1*ml1*ml1*(x2*x2)*(V2*(x5 - x6) + 2*VA*(-x2 + x5 + x6)) - 2*ml1*ml2*V2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6)))) - ml2*ml2*(V2*(-2*(x2*x2*x2*x2) + x2*x2*x2*(x3 + x4 + 4*x5) + x2*x2*(-2*x5*(x3 - x4 + 2*x5) + x1*(3*x5 - x6)) + 6*x1*x5*(x5 - x6)*x6 + 3*x1*x2*(-(x5*x5) + x6*x6)) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x3 + 2*x5)*(x5 + x6) + x2*(3*x3 + x4 + 6*x5 + 2*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)))) + ml1*ml1*(V2*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x5 + x6)*(x4 + 2*x6) + x2*(x3 + 3*x4 + 2*x5 + 6*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6))))) - 2*Enu*Mn*x1*x2*(-2*V2*x1*(x2*x2*x2*x2) + 4*V2*(x2*x2*x2*x2*x2) + V2*x1*(x2*x2*x2)*x3 + 2*VA*x1*(x2*x2*x2)*x3 - 6*V2*(x2*x2*x2*x2)*x3 - 4*VA*(x2*x2*x2*x2)*x3 + 2*V2*(x2*x2*x2)*(x3*x3) + 4*VA*(x2*x2*x2)*(x3*x3) + V2*x1*(x2*x2*x2)*x4 - 2*VA*x1*(x2*x2*x2)*x4 - 6*V2*(x2*x2*x2*x2)*x4 + 4*VA*(x2*x2*x2*x2)*x4 + 4*V2*(x2*x2*x2)*x3*x4 + 2*V2*(x2*x2*x2)*(x4*x4) - 4*VA*(x2*x2*x2)*(x4*x4) + 4*V2*(x1*x1)*(x2*x2)*x5 + 8*VA*(x1*x1)*(x2*x2)*x5 - 6*V2*x1*(x2*x2*x2)*x5 - 12*VA*x1*(x2*x2*x2)*x5 - 8*V2*(x2*x2*x2*x2)*x5 - 8*VA*(x2*x2*x2*x2)*x5 + 8*V2*x1*(x2*x2)*x3*x5 + 16*VA*x1*(x2*x2)*x3*x5 + 6*V2*(x2*x2*x2)*x3*x5 + 12*VA*(x2*x2*x2)*x3*x5 + 8*V2*x1*(x2*x2)*x4*x5 + 16*VA*x1*(x2*x2)*x4*x5 + 6*V2*(x2*x2*x2)*x4*x5 - 4*VA*(x2*x2*x2)*x4*x5 + 2*V2*(x2*x2)*x3*x4*x5 + 4*VA*(x2*x2)*x3*x4*x5 - 2*V2*(x2*x2)*(x4*x4)*x5 + 4*VA*(x2*x2)*(x4*x4)*x5 - 3*V2*(x1*x1)*x2*(x5*x5) - 6*VA*(x1*x1)*x2*(x5*x5) + 12*V2*x1*(x2*x2)*(x5*x5) + 24*VA*x1*(x2*x2)*(x5*x5) + 8*V2*(x2*x2*x2)*(x5*x5) + 16*VA*(x2*x2*x2)*(x5*x5) - 6*V2*x1*x2*x3*(x5*x5) - 12*VA*x1*x2*x3*(x5*x5) - 4*V2*(x2*x2)*x3*(x5*x5) - 8*VA*(x2*x2)*x3*(x5*x5) - 12*V2*x1*x2*x4*(x5*x5) - 24*VA*x1*x2*x4*(x5*x5) - 6*V2*x1*x2*(x5*x5*x5) - 12*VA*x1*x2*(x5*x5*x5) - 4*V2*(x2*x2)*(x5*x5*x5) - 8*VA*(x2*x2)*(x5*x5*x5) + 6*V2*x1*x4*(x5*x5*x5) + 12*VA*x1*x4*(x5*x5*x5) + 4*V2*(x1*x1)*(x2*x2)*x6 - 8*VA*(x1*x1)*(x2*x2)*x6 - 6*V2*x1*(x2*x2*x2)*x6 + 12*VA*x1*(x2*x2*x2)*x6 - 8*V2*(x2*x2*x2*x2)*x6 + 8*VA*(x2*x2*x2*x2)*x6 + 8*V2*x1*(x2*x2)*x3*x6 - 16*VA*x1*(x2*x2)*x3*x6 + 6*V2*(x2*x2*x2)*x3*x6 + 4*VA*(x2*x2*x2)*x3*x6 - 2*V2*(x2*x2)*(x3*x3)*x6 - 4*VA*(x2*x2)*(x3*x3)*x6 + 8*V2*x1*(x2*x2)*x4*x6 - 16*VA*x1*(x2*x2)*x4*x6 + 6*V2*(x2*x2*x2)*x4*x6 - 12*VA*(x2*x2*x2)*x4*x6 + 2*V2*(x2*x2)*x3*x4*x6 - 4*VA*(x2*x2)*x3*x4*x6 + 16*V2*x1*(x2*x2)*x5*x6 + 8*V2*(x2*x2*x2)*x5*x6 - 6*V2*x1*x2*x3*x5*x6 - 12*VA*x1*x2*x3*x5*x6 - 4*V2*(x2*x2)*x3*x5*x6 - 8*VA*(x2*x2)*x3*x5*x6 - 6*V2*x1*x2*x4*x5*x6 + 12*VA*x1*x2*x4*x5*x6 - 4*V2*(x2*x2)*x4*x5*x6 + 8*VA*(x2*x2)*x4*x5*x6 - 18*V2*x1*x2*(x5*x5)*x6 - 36*VA*x1*x2*(x5*x5)*x6 - 4*V2*(x2*x2)*(x5*x5)*x6 - 8*VA*(x2*x2)*(x5*x5)*x6 + 6*V2*x1*x3*(x5*x5)*x6 + 12*VA*x1*x3*(x5*x5)*x6 + 12*V2*x1*(x5*x5*x5)*x6 + 24*VA*x1*(x5*x5*x5)*x6 - 3*V2*(x1*x1)*x2*(x6*x6) + 6*VA*(x1*x1)*x2*(x6*x6) + 12*V2*x1*(x2*x2)*(x6*x6) - 24*VA*x1*(x2*x2)*(x6*x6) + 8*V2*(x2*x2*x2)*(x6*x6) - 16*VA*(x2*x2*x2)*(x6*x6) - 12*V2*x1*x2*x3*(x6*x6) + 24*VA*x1*x2*x3*(x6*x6) - 6*V2*x1*x2*x4*(x6*x6) + 12*VA*x1*x2*x4*(x6*x6) - 4*V2*(x2*x2)*x4*(x6*x6) + 8*VA*(x2*x2)*x4*(x6*x6) - 18*V2*x1*x2*x5*(x6*x6) + 36*VA*x1*x2*x5*(x6*x6) - 4*V2*(x2*x2)*x5*(x6*x6) + 8*VA*(x2*x2)*x5*(x6*x6) + 6*V2*x1*x4*x5*(x6*x6) - 12*VA*x1*x4*x5*(x6*x6) - 6*V2*x1*x2*(x6*x6*x6) + 12*VA*x1*x2*(x6*x6*x6) - 4*V2*(x2*x2)*(x6*x6*x6) + 8*VA*(x2*x2)*(x6*x6*x6) + 6*V2*x1*x3*(x6*x6*x6) - 12*VA*x1*x3*(x6*x6*x6) + 12*V2*x1*x5*(x6*x6*x6) - 24*VA*x1*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1)*ml2*V2*(x2*x2)*(-x2 + x5 + x6) + ml2*ml2*ml2*ml2*(x2*x2)*(2*VA*(x2 - x5 - x6) + V2*(-x5 + x6)) + ml1*ml1*ml1*ml1*(x2*x2)*(V2*(x5 - x6) + 2*VA*(-x2 + x5 + x6)) - 2*ml1*ml2*V2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6)))) - ml2*ml2*(V2*(-2*(x2*x2*x2*x2) + x2*x2*x2*(x3 + x4 + 4*x5) + x2*x2*(-2*x5*(x3 - x4 + 2*x5) + x1*(3*x5 - x6)) + 6*x1*x5*(x5 - x6)*x6 + 3*x1*x2*(-(x5*x5) + x6*x6)) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x3 + 2*x5)*(x5 + x6) + x2*(3*x3 + x4 + 6*x5 + 2*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)))) + ml1*ml1*(V2*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x5 + x6)*(x4 + 2*x6) + x2*(x3 + 3*x4 + 2*x5 + 6*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6))))) + x2*x2*(4*(Mn*Mn)*V2*x1*(x2*x2*x2*x2) - 8*(Mn*Mn)*V2*(x2*x2*x2*x2*x2) - 2*(Mn*Mn)*V2*x1*(x2*x2*x2)*x3 - 4*(Mn*Mn)*VA*x1*(x2*x2*x2)*x3 + 12*(Mn*Mn)*V2*(x2*x2*x2*x2)*x3 + 8*(Mn*Mn)*VA*(x2*x2*x2*x2)*x3 - 4*(Mn*Mn)*V2*(x2*x2*x2)*(x3*x3) - 8*(Mn*Mn)*VA*(x2*x2*x2)*(x3*x3) - 2*(Mn*Mn)*V2*x1*(x2*x2*x2)*x4 + 4*(Mn*Mn)*VA*x1*(x2*x2*x2)*x4 + 12*(Mn*Mn)*V2*(x2*x2*x2*x2)*x4 - 8*(Mn*Mn)*VA*(x2*x2*x2*x2)*x4 - 8*(Mn*Mn)*V2*(x2*x2*x2)*x3*x4 - 4*(Mn*Mn)*V2*(x2*x2*x2)*(x4*x4) + 8*(Mn*Mn)*VA*(x2*x2*x2)*(x4*x4) - 4*(Mn*Mn)*V2*(x1*x1)*(x2*x2)*x5 - 8*(Mn*Mn)*VA*(x1*x1)*(x2*x2)*x5 + V2*(x1*x1*x1)*(x2*x2)*x5 + 2*VA*(x1*x1*x1)*(x2*x2)*x5 + 4*(Mn*Mn)*V2*x1*(x2*x2*x2)*x5 + 8*(Mn*Mn)*VA*x1*(x2*x2*x2)*x5 - 2*V2*(x1*x1)*(x2*x2*x2)*x5 - 4*VA*(x1*x1)*(x2*x2*x2)*x5 + 16*(Mn*Mn)*V2*(x2*x2*x2*x2)*x5 + 16*(Mn*Mn)*VA*(x2*x2*x2*x2)*x5 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x3*x5 - 16*(Mn*Mn)*VA*x1*(x2*x2)*x3*x5 + 2*V2*(x1*x1)*(x2*x2)*x3*x5 + 4*VA*(x1*x1)*(x2*x2)*x3*x5 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x3*x5 - 24*(Mn*Mn)*VA*(x2*x2*x2)*x3*x5 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x4*x5 - 16*(Mn*Mn)*VA*x1*(x2*x2)*x4*x5 + 2*V2*(x1*x1)*(x2*x2)*x4*x5 + 4*VA*(x1*x1)*(x2*x2)*x4*x5 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x4*x5 + 8*(Mn*Mn)*VA*(x2*x2*x2)*x4*x5 - 4*(Mn*Mn)*V2*(x2*x2)*x3*x4*x5 - 8*(Mn*Mn)*VA*(x2*x2)*x3*x4*x5 + 4*(Mn*Mn)*V2*(x2*x2)*(x4*x4)*x5 - 8*(Mn*Mn)*VA*(x2*x2)*(x4*x4)*x5 + 2*(Mn*Mn)*V2*(x1*x1)*x2*(x5*x5) + 4*(Mn*Mn)*VA*(x1*x1)*x2*(x5*x5) - V2*(x1*x1*x1)*x2*(x5*x5) - 2*VA*(x1*x1*x1)*x2*(x5*x5) - 8*(Mn*Mn)*V2*x1*(x2*x2)*(x5*x5) - 16*(Mn*Mn)*VA*x1*(x2*x2)*(x5*x5) + 4*V2*(x1*x1)*(x2*x2)*(x5*x5) + 8*VA*(x1*x1)*(x2*x2)*(x5*x5) - 16*(Mn*Mn)*V2*(x2*x2*x2)*(x5*x5) - 32*(Mn*Mn)*VA*(x2*x2*x2)*(x5*x5) + 4*(Mn*Mn)*V2*x1*x2*x3*(x5*x5) + 8*(Mn*Mn)*VA*x1*x2*x3*(x5*x5) - 2*V2*(x1*x1)*x2*x3*(x5*x5) - 4*VA*(x1*x1)*x2*x3*(x5*x5) + 8*(Mn*Mn)*V2*(x2*x2)*x3*(x5*x5) + 16*(Mn*Mn)*VA*(x2*x2)*x3*(x5*x5) + 8*(Mn*Mn)*V2*x1*x2*x4*(x5*x5) + 16*(Mn*Mn)*VA*x1*x2*x4*(x5*x5) - 4*V2*(x1*x1)*x2*x4*(x5*x5) - 8*VA*(x1*x1)*x2*x4*(x5*x5) + 4*(Mn*Mn)*V2*x1*x2*(x5*x5*x5) + 8*(Mn*Mn)*VA*x1*x2*(x5*x5*x5) - 2*V2*(x1*x1)*x2*(x5*x5*x5) - 4*VA*(x1*x1)*x2*(x5*x5*x5) + 8*(Mn*Mn)*V2*(x2*x2)*(x5*x5*x5) + 16*(Mn*Mn)*VA*(x2*x2)*(x5*x5*x5) - 4*(Mn*Mn)*V2*x1*x4*(x5*x5*x5) - 8*(Mn*Mn)*VA*x1*x4*(x5*x5*x5) + 2*V2*(x1*x1)*x4*(x5*x5*x5) + 4*VA*(x1*x1)*x4*(x5*x5*x5) + 4*(ml1*ml1*ml1)*ml2*(Mn*Mn)*V2*(x2*x2)*(x2 - x5 - x6) - 4*(Mn*Mn)*V2*(x1*x1)*(x2*x2)*x6 + 8*(Mn*Mn)*VA*(x1*x1)*(x2*x2)*x6 + V2*(x1*x1*x1)*(x2*x2)*x6 - 2*VA*(x1*x1*x1)*(x2*x2)*x6 + 4*(Mn*Mn)*V2*x1*(x2*x2*x2)*x6 - 8*(Mn*Mn)*VA*x1*(x2*x2*x2)*x6 - 2*V2*(x1*x1)*(x2*x2*x2)*x6 + 4*VA*(x1*x1)*(x2*x2*x2)*x6 + 16*(Mn*Mn)*V2*(x2*x2*x2*x2)*x6 - 16*(Mn*Mn)*VA*(x2*x2*x2*x2)*x6 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x3*x6 + 16*(Mn*Mn)*VA*x1*(x2*x2)*x3*x6 + 2*V2*(x1*x1)*(x2*x2)*x3*x6 - 4*VA*(x1*x1)*(x2*x2)*x3*x6 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x3*x6 - 8*(Mn*Mn)*VA*(x2*x2*x2)*x3*x6 + 4*(Mn*Mn)*V2*(x2*x2)*(x3*x3)*x6 + 8*(Mn*Mn)*VA*(x2*x2)*(x3*x3)*x6 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x4*x6 + 16*(Mn*Mn)*VA*x1*(x2*x2)*x4*x6 + 2*V2*(x1*x1)*(x2*x2)*x4*x6 - 4*VA*(x1*x1)*(x2*x2)*x4*x6 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x4*x6 + 24*(Mn*Mn)*VA*(x2*x2*x2)*x4*x6 - 4*(Mn*Mn)*V2*(x2*x2)*x3*x4*x6 + 8*(Mn*Mn)*VA*(x2*x2)*x3*x4*x6 - 16*(Mn*Mn)*V2*x1*(x2*x2)*x5*x6 + 4*V2*(x1*x1)*(x2*x2)*x5*x6 - 16*(Mn*Mn)*V2*(x2*x2*x2)*x5*x6 + 4*(Mn*Mn)*V2*x1*x2*x3*x5*x6 + 8*(Mn*Mn)*VA*x1*x2*x3*x5*x6 - 2*V2*(x1*x1)*x2*x3*x5*x6 - 4*VA*(x1*x1)*x2*x3*x5*x6 + 8*(Mn*Mn)*V2*(x2*x2)*x3*x5*x6 + 16*(Mn*Mn)*VA*(x2*x2)*x3*x5*x6 + 4*(Mn*Mn)*V2*x1*x2*x4*x5*x6 - 8*(Mn*Mn)*VA*x1*x2*x4*x5*x6 - 2*V2*(x1*x1)*x2*x4*x5*x6 + 4*VA*(x1*x1)*x2*x4*x5*x6 + 8*(Mn*Mn)*V2*(x2*x2)*x4*x5*x6 - 16*(Mn*Mn)*VA*(x2*x2)*x4*x5*x6 + 12*(Mn*Mn)*V2*x1*x2*(x5*x5)*x6 + 24*(Mn*Mn)*VA*x1*x2*(x5*x5)*x6 - 6*V2*(x1*x1)*x2*(x5*x5)*x6 - 12*VA*(x1*x1)*x2*(x5*x5)*x6 + 8*(Mn*Mn)*V2*(x2*x2)*(x5*x5)*x6 + 16*(Mn*Mn)*VA*(x2*x2)*(x5*x5)*x6 - 4*(Mn*Mn)*V2*x1*x3*(x5*x5)*x6 - 8*(Mn*Mn)*VA*x1*x3*(x5*x5)*x6 + 2*V2*(x1*x1)*x3*(x5*x5)*x6 + 4*VA*(x1*x1)*x3*(x5*x5)*x6 - 8*(Mn*Mn)*V2*x1*(x5*x5*x5)*x6 - 16*(Mn*Mn)*VA*x1*(x5*x5*x5)*x6 + 4*V2*(x1*x1)*(x5*x5*x5)*x6 + 8*VA*(x1*x1)*(x5*x5*x5)*x6 + 2*(Mn*Mn)*V2*(x1*x1)*x2*(x6*x6) - 4*(Mn*Mn)*VA*(x1*x1)*x2*(x6*x6) - V2*(x1*x1*x1)*x2*(x6*x6) + 2*VA*(x1*x1*x1)*x2*(x6*x6) - 8*(Mn*Mn)*V2*x1*(x2*x2)*(x6*x6) + 16*(Mn*Mn)*VA*x1*(x2*x2)*(x6*x6) + 4*V2*(x1*x1)*(x2*x2)*(x6*x6) - 8*VA*(x1*x1)*(x2*x2)*(x6*x6) - 16*(Mn*Mn)*V2*(x2*x2*x2)*(x6*x6) + 32*(Mn*Mn)*VA*(x2*x2*x2)*(x6*x6) + 8*(Mn*Mn)*V2*x1*x2*x3*(x6*x6) - 16*(Mn*Mn)*VA*x1*x2*x3*(x6*x6) - 4*V2*(x1*x1)*x2*x3*(x6*x6) + 8*VA*(x1*x1)*x2*x3*(x6*x6) + 4*(Mn*Mn)*V2*x1*x2*x4*(x6*x6) - 8*(Mn*Mn)*VA*x1*x2*x4*(x6*x6) - 2*V2*(x1*x1)*x2*x4*(x6*x6) + 4*VA*(x1*x1)*x2*x4*(x6*x6) + 8*(Mn*Mn)*V2*(x2*x2)*x4*(x6*x6) - 16*(Mn*Mn)*VA*(x2*x2)*x4*(x6*x6) + 12*(Mn*Mn)*V2*x1*x2*x5*(x6*x6) - 24*(Mn*Mn)*VA*x1*x2*x5*(x6*x6) - 6*V2*(x1*x1)*x2*x5*(x6*x6) + 12*VA*(x1*x1)*x2*x5*(x6*x6) + 8*(Mn*Mn)*V2*(x2*x2)*x5*(x6*x6) - 16*(Mn*Mn)*VA*(x2*x2)*x5*(x6*x6) - 4*(Mn*Mn)*V2*x1*x4*x5*(x6*x6) + 8*(Mn*Mn)*VA*x1*x4*x5*(x6*x6) + 2*V2*(x1*x1)*x4*x5*(x6*x6) - 4*VA*(x1*x1)*x4*x5*(x6*x6) + 4*(Mn*Mn)*V2*x1*x2*(x6*x6*x6) - 8*(Mn*Mn)*VA*x1*x2*(x6*x6*x6) - 2*V2*(x1*x1)*x2*(x6*x6*x6) + 4*VA*(x1*x1)*x2*(x6*x6*x6) + 8*(Mn*Mn)*V2*(x2*x2)*(x6*x6*x6) - 16*(Mn*Mn)*VA*(x2*x2)*(x6*x6*x6) - 4*(Mn*Mn)*V2*x1*x3*(x6*x6*x6) + 8*(Mn*Mn)*VA*x1*x3*(x6*x6*x6) + 2*V2*(x1*x1)*x3*(x6*x6*x6) - 4*VA*(x1*x1)*x3*(x6*x6*x6) - 8*(Mn*Mn)*V2*x1*x5*(x6*x6*x6) + 16*(Mn*Mn)*VA*x1*x5*(x6*x6*x6) + 4*V2*(x1*x1)*x5*(x6*x6*x6) - 8*VA*(x1*x1)*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1*ml1)*(Mn*Mn)*(x2*x2)*(2*VA*(x2 - x5 - x6) + V2*(-x5 + x6)) + 2*(ml2*ml2*ml2*ml2)*(Mn*Mn)*(x2*x2)*(V2*(x5 - x6) + 2*VA*(-x2 + x5 + x6)) + 2*ml1*ml2*V2*(2*(ml2*ml2)*(Mn*Mn)*(x2*x2)*(x2 - x5 - x6) - x1*x1*(x2*x2*x2 - 2*(x2*x2)*(x5 + x6) - 2*x5*x6*(x5 + x6) + x2*(x5*x5 + 4*x5*x6 + x6*x6)) + 2*(Mn*Mn)*(-(x1*(x2*x2*x2 + 2*(x2*x2)*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6))))) + ml2*ml2*(x1*x1*(-(V2*(x5 - x6)*(x2*x2 + 2*x5*x6 - x2*(x5 + x6))) - 2*VA*(x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) + 2*(Mn*Mn)*(V2*(-2*(x2*x2*x2*x2) + x2*x2*x2*(x3 + x4 + 4*x5) + 2*x1*x5*(x5 - x6)*x6 + x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(-2*x5*(x3 - x4 + 2*x5) + x1*(x5 + x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x3 + 2*x5)*(x5 + x6) + x2*(3*x3 + x4 + 6*x5 + 2*x6)) + x1*(x2*x2*x2 + x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))))) + ml1*ml1*(x1*x1*(V2*(x5 - x6)*(x2*x2 + 2*x5*x6 - x2*(x5 + x6)) + 2*VA*(x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) - 2*(Mn*Mn)*(V2*(2*(x2*x2*x2*x2) + 2*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + x1*x2*(-(x5*x5) + x6*x6) - x2*x2*(2*(x3 - x4 - 2*x6)*x6 + x1*(x5 + x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x5 + x6)*(x4 + 2*x6) + x2*(x3 + 3*x4 + 2*x5 + 6*x6)) + x1*(x2*x2*x2 + x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6)))))) + A2*(2*(Enu*Enu)*(Mn*Mn)*x1*(-2*x1*(x2*x2*x2*x2) + 4*(x2*x2*x2*x2*x2) + x1*(x2*x2*x2)*x3 - 6*(x2*x2*x2*x2)*x3 + 2*(x2*x2*x2)*(x3*x3) + x1*(x2*x2*x2)*x4 - 6*(x2*x2*x2*x2)*x4 + 4*(x2*x2*x2)*x3*x4 + 2*(x2*x2*x2)*(x4*x4) + 4*(x1*x1)*(x2*x2)*x5 - 6*x1*(x2*x2*x2)*x5 - 8*(x2*x2*x2*x2)*x5 + 8*x1*(x2*x2)*x3*x5 + 6*(x2*x2*x2)*x3*x5 + 8*x1*(x2*x2)*x4*x5 + 6*(x2*x2*x2)*x4*x5 + 2*(x2*x2)*x3*x4*x5 - 2*(x2*x2)*(x4*x4)*x5 - 3*(x1*x1)*x2*(x5*x5) + 12*x1*(x2*x2)*(x5*x5) + 8*(x2*x2*x2)*(x5*x5) - 6*x1*x2*x3*(x5*x5) - 4*(x2*x2)*x3*(x5*x5) - 12*x1*x2*x4*(x5*x5) - 6*x1*x2*(x5*x5*x5) - 4*(x2*x2)*(x5*x5*x5) + 6*x1*x4*(x5*x5*x5) + 2*(ml1*ml1*ml1)*ml2*(x2*x2)*(x2 - x5 - x6) + ml1*ml1*ml1*ml1*(x2*x2)*(x5 - x6) + 4*(x1*x1)*(x2*x2)*x6 - 6*x1*(x2*x2*x2)*x6 - 8*(x2*x2*x2*x2)*x6 + 8*x1*(x2*x2)*x3*x6 + 6*(x2*x2*x2)*x3*x6 - 2*(x2*x2)*(x3*x3)*x6 + 8*x1*(x2*x2)*x4*x6 + 6*(x2*x2*x2)*x4*x6 + 2*(x2*x2)*x3*x4*x6 + 16*x1*(x2*x2)*x5*x6 + 8*(x2*x2*x2)*x5*x6 - 6*x1*x2*x3*x5*x6 - 4*(x2*x2)*x3*x5*x6 - 6*x1*x2*x4*x5*x6 - 4*(x2*x2)*x4*x5*x6 - 18*x1*x2*(x5*x5)*x6 - 4*(x2*x2)*(x5*x5)*x6 + 6*x1*x3*(x5*x5)*x6 + 12*x1*(x5*x5*x5)*x6 - 3*(x1*x1)*x2*(x6*x6) + 12*x1*(x2*x2)*(x6*x6) + 8*(x2*x2*x2)*(x6*x6) - 12*x1*x2*x3*(x6*x6) - 6*x1*x2*x4*(x6*x6) - 4*(x2*x2)*x4*(x6*x6) - 18*x1*x2*x5*(x6*x6) - 4*(x2*x2)*x5*(x6*x6) + 6*x1*x4*x5*(x6*x6) - 6*x1*x2*(x6*x6*x6) - 4*(x2*x2)*(x6*x6*x6) + 6*x1*x3*(x6*x6*x6) + 12*x1*x5*(x6*x6*x6) + ml2*ml2*ml2*ml2*(x2*x2)*(-x5 + x6) + ml2*ml2*(2*(x2*x2*x2*x2) - x2*x2*x2*(x3 + x4 + 4*x5) + 6*x1*x5*x6*(-x5 + x6) + 3*x1*x2*(x5*x5 - x6*x6) + x2*x2*(2*x5*(x3 - x4 + 2*x5) + x1*(-3*x5 + x6))) + ml1*ml1*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*ml1*ml2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6))))) - 2*Enu*Mn*x1*x2*(-2*x1*(x2*x2*x2*x2) + 4*(x2*x2*x2*x2*x2) + x1*(x2*x2*x2)*x3 - 6*(x2*x2*x2*x2)*x3 + 2*(x2*x2*x2)*(x3*x3) + x1*(x2*x2*x2)*x4 - 6*(x2*x2*x2*x2)*x4 + 4*(x2*x2*x2)*x3*x4 + 2*(x2*x2*x2)*(x4*x4) + 4*(x1*x1)*(x2*x2)*x5 - 6*x1*(x2*x2*x2)*x5 - 8*(x2*x2*x2*x2)*x5 + 8*x1*(x2*x2)*x3*x5 + 6*(x2*x2*x2)*x3*x5 + 8*x1*(x2*x2)*x4*x5 + 6*(x2*x2*x2)*x4*x5 + 2*(x2*x2)*x3*x4*x5 - 2*(x2*x2)*(x4*x4)*x5 - 3*(x1*x1)*x2*(x5*x5) + 12*x1*(x2*x2)*(x5*x5) + 8*(x2*x2*x2)*(x5*x5) - 6*x1*x2*x3*(x5*x5) - 4*(x2*x2)*x3*(x5*x5) - 12*x1*x2*x4*(x5*x5) - 6*x1*x2*(x5*x5*x5) - 4*(x2*x2)*(x5*x5*x5) + 6*x1*x4*(x5*x5*x5) + 2*(ml1*ml1*ml1)*ml2*(x2*x2)*(x2 - x5 - x6) + ml1*ml1*ml1*ml1*(x2*x2)*(x5 - x6) + 4*(x1*x1)*(x2*x2)*x6 - 6*x1*(x2*x2*x2)*x6 - 8*(x2*x2*x2*x2)*x6 + 8*x1*(x2*x2)*x3*x6 + 6*(x2*x2*x2)*x3*x6 - 2*(x2*x2)*(x3*x3)*x6 + 8*x1*(x2*x2)*x4*x6 + 6*(x2*x2*x2)*x4*x6 + 2*(x2*x2)*x3*x4*x6 + 16*x1*(x2*x2)*x5*x6 + 8*(x2*x2*x2)*x5*x6 - 6*x1*x2*x3*x5*x6 - 4*(x2*x2)*x3*x5*x6 - 6*x1*x2*x4*x5*x6 - 4*(x2*x2)*x4*x5*x6 - 18*x1*x2*(x5*x5)*x6 - 4*(x2*x2)*(x5*x5)*x6 + 6*x1*x3*(x5*x5)*x6 + 12*x1*(x5*x5*x5)*x6 - 3*(x1*x1)*x2*(x6*x6) + 12*x1*(x2*x2)*(x6*x6) + 8*(x2*x2*x2)*(x6*x6) - 12*x1*x2*x3*(x6*x6) - 6*x1*x2*x4*(x6*x6) - 4*(x2*x2)*x4*(x6*x6) - 18*x1*x2*x5*(x6*x6) - 4*(x2*x2)*x5*(x6*x6) + 6*x1*x4*x5*(x6*x6) - 6*x1*x2*(x6*x6*x6) - 4*(x2*x2)*(x6*x6*x6) + 6*x1*x3*(x6*x6*x6) + 12*x1*x5*(x6*x6*x6) + ml2*ml2*ml2*ml2*(x2*x2)*(-x5 + x6) + ml2*ml2*(2*(x2*x2*x2*x2) - x2*x2*x2*(x3 + x4 + 4*x5) + 6*x1*x5*x6*(-x5 + x6) + 3*x1*x2*(x5*x5 - x6*x6) + x2*x2*(2*x5*(x3 - x4 + 2*x5) + x1*(-3*x5 + x6))) + ml1*ml1*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*ml1*ml2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6))))) + x2*x2*(-4*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2*x2) + 4*(Mn*Mn)*x1*(x2*x2*x2*x2) - 8*(Mn*Mn)*(x2*x2*x2*x2*x2) + 2*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2)*x3 - 2*(Mn*Mn)*x1*(x2*x2*x2)*x3 + 12*(Mn*Mn)*(x2*x2*x2*x2)*x3 - 4*(Mn*Mn)*(x2*x2*x2)*(x3*x3) + 2*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2)*x4 - 2*(Mn*Mn)*x1*(x2*x2*x2)*x4 + 12*(Mn*Mn)*(x2*x2*x2*x2)*x4 - 8*(Mn*Mn)*(x2*x2*x2)*x3*x4 - 4*(Mn*Mn)*(x2*x2*x2)*(x4*x4) + 2*(ml2*ml2*ml2*ml2)*(Mn*Mn)*(x2*x2)*x5 + 2*(ml2*ml2)*(Mn*Mn)*x1*(x2*x2)*x5 - ml2*ml2*(x1*x1)*(x2*x2)*x5 - 4*(Mn*Mn)*(x1*x1)*(x2*x2)*x5 + x1*x1*x1*(x2*x2)*x5 + 8*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2)*x5 + 4*(Mn*Mn)*x1*(x2*x2*x2)*x5 - 2*(x1*x1)*(x2*x2*x2)*x5 + 16*(Mn*Mn)*(x2*x2*x2*x2)*x5 - 4*(ml2*ml2)*(Mn*Mn)*(x2*x2)*x3*x5 - 8*(Mn*Mn)*x1*(x2*x2)*x3*x5 + 2*(x1*x1)*(x2*x2)*x3*x5 - 12*(Mn*Mn)*(x2*x2*x2)*x3*x5 + 4*(ml2*ml2)*(Mn*Mn)*(x2*x2)*x4*x5 - 8*(Mn*Mn)*x1*(x2*x2)*x4*x5 + 2*(x1*x1)*(x2*x2)*x4*x5 - 12*(Mn*Mn)*(x2*x2*x2)*x4*x5 - 4*(Mn*Mn)*(x2*x2)*x3*x4*x5 + 4*(Mn*Mn)*(x2*x2)*(x4*x4)*x5 - 2*(ml2*ml2)*(Mn*Mn)*x1*x2*(x5*x5) + ml2*ml2*(x1*x1)*x2*(x5*x5) + 2*(Mn*Mn)*(x1*x1)*x2*(x5*x5) - x1*x1*x1*x2*(x5*x5) - 8*(ml2*ml2)*(Mn*Mn)*(x2*x2)*(x5*x5) - 8*(Mn*Mn)*x1*(x2*x2)*(x5*x5) + 4*(x1*x1)*(x2*x2)*(x5*x5) - 16*(Mn*Mn)*(x2*x2*x2)*(x5*x5) + 4*(Mn*Mn)*x1*x2*x3*(x5*x5) - 2*(x1*x1)*x2*x3*(x5*x5) + 8*(Mn*Mn)*(x2*x2)*x3*(x5*x5) + 8*(Mn*Mn)*x1*x2*x4*(x5*x5) - 4*(x1*x1)*x2*x4*(x5*x5) + 4*(Mn*Mn)*x1*x2*(x5*x5*x5) - 2*(x1*x1)*x2*(x5*x5*x5) + 8*(Mn*Mn)*(x2*x2)*(x5*x5*x5) - 4*(Mn*Mn)*x1*x4*(x5*x5*x5) + 2*(x1*x1)*x4*(x5*x5*x5) - 2*(ml2*ml2*ml2*ml2)*(Mn*Mn)*(x2*x2)*x6 + 2*(ml2*ml2)*(Mn*Mn)*x1*(x2*x2)*x6 + ml2*ml2*(x1*x1)*(x2*x2)*x6 - 4*(Mn*Mn)*(x1*x1)*(x2*x2)*x6 + x1*x1*x1*(x2*x2)*x6 + 4*(Mn*Mn)*x1*(x2*x2*x2)*x6 - 2*(x1*x1)*(x2*x2*x2)*x6 + 16*(Mn*Mn)*(x2*x2*x2*x2)*x6 - 8*(Mn*Mn)*x1*(x2*x2)*x3*x6 + 2*(x1*x1)*(x2*x2)*x3*x6 - 12*(Mn*Mn)*(x2*x2*x2)*x3*x6 + 4*(Mn*Mn)*(x2*x2)*(x3*x3)*x6 - 8*(Mn*Mn)*x1*(x2*x2)*x4*x6 + 2*(x1*x1)*(x2*x2)*x4*x6 - 12*(Mn*Mn)*(x2*x2*x2)*x4*x6 - 4*(Mn*Mn)*(x2*x2)*x3*x4*x6 - 16*(Mn*Mn)*x1*(x2*x2)*x5*x6 + 4*(x1*x1)*(x2*x2)*x5*x6 - 16*(Mn*Mn)*(x2*x2*x2)*x5*x6 + 4*(Mn*Mn)*x1*x2*x3*x5*x6 - 2*(x1*x1)*x2*x3*x5*x6 + 8*(Mn*Mn)*(x2*x2)*x3*x5*x6 + 4*(Mn*Mn)*x1*x2*x4*x5*x6 - 2*(x1*x1)*x2*x4*x5*x6 + 8*(Mn*Mn)*(x2*x2)*x4*x5*x6 + 4*(ml2*ml2)*(Mn*Mn)*x1*(x5*x5)*x6 - 2*(ml2*ml2)*(x1*x1)*(x5*x5)*x6 + 12*(Mn*Mn)*x1*x2*(x5*x5)*x6 - 6*(x1*x1)*x2*(x5*x5)*x6 + 8*(Mn*Mn)*(x2*x2)*(x5*x5)*x6 - 4*(Mn*Mn)*x1*x3*(x5*x5)*x6 + 2*(x1*x1)*x3*(x5*x5)*x6 - 8*(Mn*Mn)*x1*(x5*x5*x5)*x6 + 4*(x1*x1)*(x5*x5*x5)*x6 + 2*(ml2*ml2)*(Mn*Mn)*x1*x2*(x6*x6) - ml2*ml2*(x1*x1)*x2*(x6*x6) + 2*(Mn*Mn)*(x1*x1)*x2*(x6*x6) - x1*x1*x1*x2*(x6*x6) - 8*(Mn*Mn)*x1*(x2*x2)*(x6*x6) + 4*(x1*x1)*(x2*x2)*(x6*x6) - 16*(Mn*Mn)*(x2*x2*x2)*(x6*x6) + 8*(Mn*Mn)*x1*x2*x3*(x6*x6) - 4*(x1*x1)*x2*x3*(x6*x6) + 4*(Mn*Mn)*x1*x2*x4*(x6*x6) - 2*(x1*x1)*x2*x4*(x6*x6) + 8*(Mn*Mn)*(x2*x2)*x4*(x6*x6) - 4*(ml2*ml2)*(Mn*Mn)*x1*x5*(x6*x6) + 2*(ml2*ml2)*(x1*x1)*x5*(x6*x6) + 12*(Mn*Mn)*x1*x2*x5*(x6*x6) - 6*(x1*x1)*x2*x5*(x6*x6) + 8*(Mn*Mn)*(x2*x2)*x5*(x6*x6) - 4*(Mn*Mn)*x1*x4*x5*(x6*x6) + 2*(x1*x1)*x4*x5*(x6*x6) + 4*(Mn*Mn)*x1*x2*(x6*x6*x6) - 2*(x1*x1)*x2*(x6*x6*x6) + 8*(Mn*Mn)*(x2*x2)*(x6*x6*x6) - 4*(Mn*Mn)*x1*x3*(x6*x6*x6) + 2*(x1*x1)*x3*(x6*x6*x6) - 8*(Mn*Mn)*x1*x5*(x6*x6*x6) + 4*(x1*x1)*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1*ml1)*(Mn*Mn)*(x2*x2)*(-x5 + x6) + 4*(ml1*ml1*ml1)*ml2*(Mn*Mn)*(x2*x2)*(-x2 + x5 + x6) + ml1*ml1*(x1*x1*(x5 - x6)*(x2*x2 + 2*x5*x6 - x2*(x5 + x6)) + 2*(Mn*Mn)*(-2*(x2*x2*x2*x2) + 2*x1*x5*x6*(-x5 + x6) + x2*x2*x2*(x3 + x4 + 4*x6) + x1*x2*(x5*x5 - x6*x6) + x2*x2*(2*(x3 - x4 - 2*x6)*x6 + x1*(x5 + x6)))) - 2*ml1*ml2*(2*(ml2*ml2)*(Mn*Mn)*(x2*x2)*(x2 - x5 - x6) - x1*x1*(x2*x2*x2 - 2*(x2*x2)*(x5 + x6) - 2*x5*x6*(x5 + x6) + x2*(x5*x5 + 4*x5*x6 + x6*x6)) + 2*(Mn*Mn)*(-(x1*(x2*x2*x2 + 2*(x2*x2)*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6)))))))))*(Z*Z))/(Enu*Enu*(Mn*Mn)*(x1*x1)*(x2*x2*x2*x2)*((x1 + 2*x3)*(x1 + 2*x3))*((x1 + 2*x4)*(x1 + 2*x4)));

  ///////////////////////////////////////////////////////////////////////////////////
  // JACOBIANS AND PREFACTORS 
  dsigma *= Jacob;

  // Changing units to zeptobarn = zb = 1e-45 cm^2
  dsigma *= (GeV2_to_cm2*1e45);
  

  myMC->xphys[0] = x1;
  myMC->xphys[1] = x2;
  myMC->xphys[2] = x3;
  myMC->xphys[3] = x4;
  myMC->xphys[4] = x5;
  myMC->xphys[5] = x6;
  myMC->xphys[6] = x7;
  myMC->xphys[7] = x8;
  myMC->xphys[8] = Enu;
  myMC->xphys[9] = dsigma;


  return 0;
}

int CohFromIntToPhysical(const cubareal xx[], void *MC){

  long double u1_s    = xx[0];
  long double u2_s    = xx[1];
  long double u3_s    = xx[2];
  long double PHI2_s  = xx[3];
  long double x5_s    = xx[4];
  long double u6_s    = xx[5];
  long double x7_s    = xx[6];
  long double x8_s    = xx[7];
  long double Enu_s   = xx[8];

  // Random angles needed to specify 4 vectors
  long double x7 = x7_s*2*M_PI;
  long double x8 = x8_s*2*M_PI;

  //////////////////////////////////////////////////
  // get data from user

  tridentMC* myMC = (tridentMC *)MC;
  
  long double ml1 = myMC->ml1;
  long double ml2 = myMC->ml2;
  long double A = myMC->A;
  long double Z = myMC->Z;
  long double Mn = myMC->Mn;
  long double mzprime = myMC->mzprime;

  long double Diag11 = myMC->terms[0];
  long double Diag22 = myMC->terms[1];
  long double Diag12 = myMC->terms[2];
  long double V2     = myMC->terms[3];
  long double A2     = myMC->terms[4];
  long double VA     = myMC->terms[5];
  long double BSM    = myMC->terms[6];

  long double Enu = Enu_s*(myMC->Emax - myMC->Emin) + myMC->Emin;


  //////////////////////////////////
  long double x1_u = 2*SQR(Enu) / ( 1 + 2 * Enu/Mn ) * ( 1 - SQR(ml1+ml2)/2.0/SQR(Enu) * (1 + Enu/Mn) +
                  sqrt( SQR(1 - SQR(ml1 + ml2)/2.0/SQR(Enu)*(1+Enu/Mn) ) 
                    - pow(ml1 + ml2, 4)/4.0/pow(Enu,4)*(1 + 2*Enu/Mn) )  );
  long double x1_l = pow(ml1+ml2,4)/(x1_u)/(1.0+2*Enu/Mn);
  
  long double u1 = u1_s*(log(x1_u)-log(x1_l)) + log(x1_l);

  long double x1 = exp(u1);
  // std::cout<<x1_l<<std::endl;

  //////////////////////////////////
  long double x2_l = 1.0/2.0 * ( x1 + SQR(ml1 + ml2) );
  long double x2_u = Enu * (sqrt( x1 + SQR(x1)/4.0/SQR(Mn) ) - x1/2.0/Mn );
  // long double x2 = x2_s*(x2_u-x2_l) + x2_l;

  long double u2_u = (1.0+2.0*Enu/Mn)*(x1_u-x1)*(x1-x1_l)  / ( SQR(ml1 + ml2) + x1*(1+Enu/Mn) + 2*Enu*sqrt(x1+SQR(x1)/4.0/SQR(Mn)) );

  long double u2 = u2_s*(u2_u-0.0) + 0.0;

  long double x2 = 0.5*(u2 + SQR(ml1 + ml2) + x1);

  /////////////////
  long double x5_l = x2/2.0/(2*x2-x1) * (  2*x2 - x1 + SQR(ml2) - SQR(ml1)  
        - sqrt( SQR( 2*x2-x1 + SQR(ml2) - SQR(ml1) ) - 4*SQR(ml2)*(2*x2-x1) ) );
  
  long double x5_u = x2/2.0/(2*x2-x1) * (  2*x2 - x1 + SQR(ml2) - SQR(ml1)  
        + sqrt( SQR( 2*x2-x1 + SQR(ml2) - SQR(ml1) ) - 4*SQR(ml2)*(2*x2-x1) ) );
  

  long double x5 = x5_s*(x5_u-x5_l) + x5_l; 

  ////////////////
  long double x3_l = SQR(ml2)/2.0 * x2/x5 - x1/2.0 * x5/x2;

  long double x3_u = 0.5*(2*x2 - x1 + SQR(ml2) - SQR(ml1) )  - x5;

  long double u3 = u3_s*(log(2*x3_u + x1)-log(2*x3_l + x1)) + log(2*x3_l + x1);

  long double x3 = 0.5*(exp(u3)- x1);


  ////////////////
  // Useful Definitions in Frame p1vec + qvec - p3vec = 0
  ////////////////
  long double Wc2 = 2*(x2-x3-x5) - x1 + SQR(ml2);

  long double E1 =  (x2 - x5)/sqrt(Wc2) ;
  long double E4 =  (Wc2+SQR(ml1))/2.0/sqrt(Wc2) ;
  long double q0 =  (x2-x1-x3)/sqrt(Wc2) ;
  long double qvec = sqrt(SQR(q0) + x1)  ;
  long double p4vec = (Wc2 - SQR(ml1))/2.0/sqrt(Wc2)  ;

  long double Cq = (q0*E1 - x2)/qvec/E1;
  long double Sq = sqrt(1 - SQR(Cq));

  long double E2 = (Wc2 - SQR(ml1))/2/sqrt(Wc2);


  ////////////////

  long double m6, u6, u6_l, u6_u, prop, jacob_u6;
  long double m6_l = 0;
  long double m6_u = 2*E1*E2;

  switch ((int) BSM) 
  {
    case (SMonly):
      u6_l = m6_l;
      u6_u = m6_u;
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6 = u6;
      jacob_u6 = 1.0;

      V2     = myMC->terms[3];
      A2     = myMC->terms[4];
      VA     = myMC->terms[5];
      break;

    case (INTERFERENCE):

      //////////////////////////////////////////////////////////
      // BSM PROPAGATOR 
      // long double prop = (1.0/(-2*(m6/mzprime/mzprime) - 1.0))/mzprime/mzprime;
      u6_l = log(1.0/(2.0*m6_u + mzprime*mzprime));
      u6_u = log(1.0/(2.0*m6_l + mzprime*mzprime));
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6   = -mzprime*mzprime/2.0 + exp(-u6)/2.0;

      prop = -exp(u6);
      jacob_u6 = exp(-u6)/2.0;

      V2     = myMC->terms[3] * prop;
      A2     = myMC->terms[4] * prop;
      VA     = myMC->terms[5] * prop;
      break;

    case (BSMonly):

      u6_l = 1.0/(2.0*m6_u + mzprime*mzprime);
      u6_u = 1.0/(2.0*m6_l + mzprime*mzprime);
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6   = -mzprime*mzprime/2.0 + 1.0/2.0/u6;

      prop = -u6;
      jacob_u6 = 1.0/(u6*u6)/2.0;

      V2     = myMC->terms[3] * prop*prop;
      A2     = myMC->terms[4] * prop*prop;
      VA     = myMC->terms[5] * prop*prop;
      break;
    
    case (SMandBSM):
     
      u6_l = m6_l;
      u6_u = m6_u;
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6 = u6;
      jacob_u6 = 1.0;

      prop = 1.0/(2.0*m6 + mzprime*mzprime);

      V2     = (myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop) ;
      A2     = (myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop) ;
      VA     = (myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop);
      break;

    default:
      std::cout << "Could not determine what contributions (BSM) to include." << std::endl;
      break;
  }


  long double C2 = 1.0 - m6 / E1 / E2;
  long double S2 = sqrt(1.0 - SQR(C2));


  ////////////////
  long double PHI2 = PHI2_s*(2*M_PI-0.0) + 0.0;

  long double x4 = q0*E4 + (Sq*S2*cos(PHI2) + Cq*C2)*E2*qvec;

  long double C4 = (E4*q0 - x4)/(p4vec*qvec);
  long double S4 = sqrt(1.0 - SQR(C4));

  long double x6 = E1*E4 + E1*E2*C2; 

  //////////////////////////////////////////////////////////
  // Multiply by the appropriate Jacobians for our invariants and phase space factors
  long double Jacob =  (Wc2 - SQR(ml1))/8.0/Wc2/E1/E2/x2/pow(2*M_PI,6.0);

  // Jacobian due to change from CSW(+Matheus) to better integration variables AND VEGAS
  Jacob *=  (log(x1_u)-log(x1_l))*exp(u1)*
            (u2_u-0.0)*0.5*
            (log(2*x3_u + x1)-log(2*x3_l + x1))*0.5*exp(u3)*
            (2*M_PI)*
            (x5_u-x5_l)*
            (u6_u-u6_l)*jacob_u6*
            (2*M_PI)*
            (2*M_PI);


  //////////////////////////////////////////////////////////
  // MATRIX ELEMENT ITSELF
  long double FormFactor = FF_WS(sqrt(x1), A);
  // includes Z*Z for the coherent enhancement
  long double dsigma = (8*(alphaQED*alphaQED)*(FormFactor*FormFactor)*(Gf*Gf)*(Diag22*((x1 + 2*x4)*(x1 + 2*x4))*(A2*(-2*(ml2*ml2)*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) - 2*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) + 2*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5)*x5 - x1*x1*(x2*x2)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 3*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 4*ml1*(ml2*ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(x2 - x5 - x6) + 2*ml1*ml2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6) - 4*ml1*ml2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(x2 - x5 - x6) - 2*ml1*ml2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*(x2 - x5 - x6) - 4*ml1*ml2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5)*(x2 - x5 - x6) - 4*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x6 - 4*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(-x1 + x2 - x3 - x4)*x6 - 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*x6 + 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 + x1*x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6)*x6 + 2*Mn*(x1*x1)*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x3*x5*(x2 - x5 - x6)*x6 + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x3*x5*(x2 - x5 - x6)*x6 - 2*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + x1*x1*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) + V2*(-2*(ml2*ml2)*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) - 2*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) + 2*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5)*x5 - x1*x1*(x2*x2)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 3*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 4*ml1*(ml2*ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(x2 - x5 - x6) - 2*ml1*ml2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6) + 4*ml1*ml2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(x2 - x5 - x6) + 2*ml1*ml2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*(x2 - x5 - x6) + 4*ml1*ml2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5)*(x2 - x5 - x6) - 4*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x6 - 4*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(-x1 + x2 - x3 - x4)*x6 - 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*x6 + 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 + x1*x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6)*x6 + 2*Mn*(x1*x1)*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x3*x5*(x2 - x5 - x6)*x6 + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x3*x5*(x2 - x5 - x6)*x6 - 2*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + x1*x1*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) + 2*VA*(-2*(ml2*ml2)*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) - 2*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) + 2*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5)*x5 - x1*x1*(x2*x2)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 3*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 4*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x6 + 4*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x3*(-x1 + x2 - x3 - x4)*x6 + 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*x6 - 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 - 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 - 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 - x1*x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6)*x6 - 2*Mn*(x1*x1)*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(x2 - x5 - x6)*x6 - 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x3*x5*(x2 - x5 - x6)*x6 - 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x3*x5*(x2 - x5 - x6)*x6 + Mn*x1*(x2*x2)*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*(ml2*ml2)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - x1*x1*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6))) + Diag11*((x1 + 2*x3)*(x1 + 2*x3))*(A2*(-4*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x5 - 4*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x4*x5 - 2*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 2*(ml1*ml1)*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) - 2*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x4*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) - 4*(ml1*ml1*ml1)*ml2*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x4*(x2 - x5 - x6) - 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*x6 + 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 + x1*x1*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 - 2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 + 2*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6)*x6 + 2*ml1*ml2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x2 - x5 - x6)*x6 - 4*ml1*ml2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6)*x6 + x1*x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6)*x6 + 2*Mn*(x1*x1)*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x4*x5*(x2 - x5 - x6)*x6 + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x4*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) - 2*ml1*ml2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x2 - x5 - x6)*(x6*x6) - 4*ml1*ml2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x2 - x5 - x6)*(x6*x6) - x1*x1*(x2*x2)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 3*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) + V2*(-4*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x5 - 4*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x4*x5 - 2*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 2*(ml1*ml1)*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) - 2*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x4*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) + 4*(ml1*ml1*ml1)*ml2*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x4*(x2 - x5 - x6) - 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*x6 + 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 + x1*x1*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 - 2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 + 2*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6)*x6 - 2*ml1*ml2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x2 - x5 - x6)*x6 + 4*ml1*ml2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x2 - x5 - x6)*x6 + x1*x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6)*x6 + 2*Mn*(x1*x1)*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x4*x5*(x2 - x5 - x6)*x6 + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x4*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) + 2*ml1*ml2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x2 - x5 - x6)*(x6*x6) + 4*ml1*ml2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x2 - x5 - x6)*(x6*x6) - x1*x1*(x2*x2)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 3*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) + 2*VA*(-4*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x5 - 4*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(-x1 + x2 - x3 - x4)*x4*x5 - 2*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*(ml1*ml1)*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) + 2*Mn*(x2*x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x4*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) - 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*x6 + 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 + x1*x1*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 - 2*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 - 2*(ml1*ml1)*Mn*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6)*x6 + x1*x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(x2 - x5 - x6)*x6 + 2*Mn*(x1*x1)*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x4*x5*(x2 - x5 - x6)*x6 + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x4*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) + 4*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) - x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) - 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) + x1*x1*(x2*x2)*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 3*Mn*x1*(x2*x2)*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*(x1*x1)*x2*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 4*Mn*x1*x2*(-(Enu*Enu*Mn*x1) + Enu*x1*x2 + Mn*(x2*x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + x1*x1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*Mn*x1*(Enu*Enu*Mn*x1 - Enu*x1*x2 - Mn*(x2*x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6))) + Diag12*(x1 + 2*x3)*(x1 + 2*x4)*(2*(Enu*Enu)*(Mn*Mn)*x1*(-2*V2*x1*(x2*x2*x2*x2) + 4*V2*(x2*x2*x2*x2*x2) + V2*x1*(x2*x2*x2)*x3 + 2*VA*x1*(x2*x2*x2)*x3 - 6*V2*(x2*x2*x2*x2)*x3 - 4*VA*(x2*x2*x2*x2)*x3 + 2*V2*(x2*x2*x2)*(x3*x3) + 4*VA*(x2*x2*x2)*(x3*x3) + V2*x1*(x2*x2*x2)*x4 - 2*VA*x1*(x2*x2*x2)*x4 - 6*V2*(x2*x2*x2*x2)*x4 + 4*VA*(x2*x2*x2*x2)*x4 + 4*V2*(x2*x2*x2)*x3*x4 + 2*V2*(x2*x2*x2)*(x4*x4) - 4*VA*(x2*x2*x2)*(x4*x4) + 4*V2*(x1*x1)*(x2*x2)*x5 + 8*VA*(x1*x1)*(x2*x2)*x5 - 6*V2*x1*(x2*x2*x2)*x5 - 12*VA*x1*(x2*x2*x2)*x5 - 8*V2*(x2*x2*x2*x2)*x5 - 8*VA*(x2*x2*x2*x2)*x5 + 8*V2*x1*(x2*x2)*x3*x5 + 16*VA*x1*(x2*x2)*x3*x5 + 6*V2*(x2*x2*x2)*x3*x5 + 12*VA*(x2*x2*x2)*x3*x5 + 8*V2*x1*(x2*x2)*x4*x5 + 16*VA*x1*(x2*x2)*x4*x5 + 6*V2*(x2*x2*x2)*x4*x5 - 4*VA*(x2*x2*x2)*x4*x5 + 2*V2*(x2*x2)*x3*x4*x5 + 4*VA*(x2*x2)*x3*x4*x5 - 2*V2*(x2*x2)*(x4*x4)*x5 + 4*VA*(x2*x2)*(x4*x4)*x5 - 3*V2*(x1*x1)*x2*(x5*x5) - 6*VA*(x1*x1)*x2*(x5*x5) + 12*V2*x1*(x2*x2)*(x5*x5) + 24*VA*x1*(x2*x2)*(x5*x5) + 8*V2*(x2*x2*x2)*(x5*x5) + 16*VA*(x2*x2*x2)*(x5*x5) - 6*V2*x1*x2*x3*(x5*x5) - 12*VA*x1*x2*x3*(x5*x5) - 4*V2*(x2*x2)*x3*(x5*x5) - 8*VA*(x2*x2)*x3*(x5*x5) - 12*V2*x1*x2*x4*(x5*x5) - 24*VA*x1*x2*x4*(x5*x5) - 6*V2*x1*x2*(x5*x5*x5) - 12*VA*x1*x2*(x5*x5*x5) - 4*V2*(x2*x2)*(x5*x5*x5) - 8*VA*(x2*x2)*(x5*x5*x5) + 6*V2*x1*x4*(x5*x5*x5) + 12*VA*x1*x4*(x5*x5*x5) + 4*V2*(x1*x1)*(x2*x2)*x6 - 8*VA*(x1*x1)*(x2*x2)*x6 - 6*V2*x1*(x2*x2*x2)*x6 + 12*VA*x1*(x2*x2*x2)*x6 - 8*V2*(x2*x2*x2*x2)*x6 + 8*VA*(x2*x2*x2*x2)*x6 + 8*V2*x1*(x2*x2)*x3*x6 - 16*VA*x1*(x2*x2)*x3*x6 + 6*V2*(x2*x2*x2)*x3*x6 + 4*VA*(x2*x2*x2)*x3*x6 - 2*V2*(x2*x2)*(x3*x3)*x6 - 4*VA*(x2*x2)*(x3*x3)*x6 + 8*V2*x1*(x2*x2)*x4*x6 - 16*VA*x1*(x2*x2)*x4*x6 + 6*V2*(x2*x2*x2)*x4*x6 - 12*VA*(x2*x2*x2)*x4*x6 + 2*V2*(x2*x2)*x3*x4*x6 - 4*VA*(x2*x2)*x3*x4*x6 + 16*V2*x1*(x2*x2)*x5*x6 + 8*V2*(x2*x2*x2)*x5*x6 - 6*V2*x1*x2*x3*x5*x6 - 12*VA*x1*x2*x3*x5*x6 - 4*V2*(x2*x2)*x3*x5*x6 - 8*VA*(x2*x2)*x3*x5*x6 - 6*V2*x1*x2*x4*x5*x6 + 12*VA*x1*x2*x4*x5*x6 - 4*V2*(x2*x2)*x4*x5*x6 + 8*VA*(x2*x2)*x4*x5*x6 - 18*V2*x1*x2*(x5*x5)*x6 - 36*VA*x1*x2*(x5*x5)*x6 - 4*V2*(x2*x2)*(x5*x5)*x6 - 8*VA*(x2*x2)*(x5*x5)*x6 + 6*V2*x1*x3*(x5*x5)*x6 + 12*VA*x1*x3*(x5*x5)*x6 + 12*V2*x1*(x5*x5*x5)*x6 + 24*VA*x1*(x5*x5*x5)*x6 - 3*V2*(x1*x1)*x2*(x6*x6) + 6*VA*(x1*x1)*x2*(x6*x6) + 12*V2*x1*(x2*x2)*(x6*x6) - 24*VA*x1*(x2*x2)*(x6*x6) + 8*V2*(x2*x2*x2)*(x6*x6) - 16*VA*(x2*x2*x2)*(x6*x6) - 12*V2*x1*x2*x3*(x6*x6) + 24*VA*x1*x2*x3*(x6*x6) - 6*V2*x1*x2*x4*(x6*x6) + 12*VA*x1*x2*x4*(x6*x6) - 4*V2*(x2*x2)*x4*(x6*x6) + 8*VA*(x2*x2)*x4*(x6*x6) - 18*V2*x1*x2*x5*(x6*x6) + 36*VA*x1*x2*x5*(x6*x6) - 4*V2*(x2*x2)*x5*(x6*x6) + 8*VA*(x2*x2)*x5*(x6*x6) + 6*V2*x1*x4*x5*(x6*x6) - 12*VA*x1*x4*x5*(x6*x6) - 6*V2*x1*x2*(x6*x6*x6) + 12*VA*x1*x2*(x6*x6*x6) - 4*V2*(x2*x2)*(x6*x6*x6) + 8*VA*(x2*x2)*(x6*x6*x6) + 6*V2*x1*x3*(x6*x6*x6) - 12*VA*x1*x3*(x6*x6*x6) + 12*V2*x1*x5*(x6*x6*x6) - 24*VA*x1*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1)*ml2*V2*(x2*x2)*(-x2 + x5 + x6) + ml2*ml2*ml2*ml2*(x2*x2)*(2*VA*(x2 - x5 - x6) + V2*(-x5 + x6)) + ml1*ml1*ml1*ml1*(x2*x2)*(V2*(x5 - x6) + 2*VA*(-x2 + x5 + x6)) - 2*ml1*ml2*V2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6)))) - ml2*ml2*(V2*(-2*(x2*x2*x2*x2) + x2*x2*x2*(x3 + x4 + 4*x5) + x2*x2*(-2*x5*(x3 - x4 + 2*x5) + x1*(3*x5 - x6)) + 6*x1*x5*(x5 - x6)*x6 + 3*x1*x2*(-(x5*x5) + x6*x6)) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x3 + 2*x5)*(x5 + x6) + x2*(3*x3 + x4 + 6*x5 + 2*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)))) + ml1*ml1*(V2*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x5 + x6)*(x4 + 2*x6) + x2*(x3 + 3*x4 + 2*x5 + 6*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6))))) - 2*Enu*Mn*x1*x2*(-2*V2*x1*(x2*x2*x2*x2) + 4*V2*(x2*x2*x2*x2*x2) + V2*x1*(x2*x2*x2)*x3 + 2*VA*x1*(x2*x2*x2)*x3 - 6*V2*(x2*x2*x2*x2)*x3 - 4*VA*(x2*x2*x2*x2)*x3 + 2*V2*(x2*x2*x2)*(x3*x3) + 4*VA*(x2*x2*x2)*(x3*x3) + V2*x1*(x2*x2*x2)*x4 - 2*VA*x1*(x2*x2*x2)*x4 - 6*V2*(x2*x2*x2*x2)*x4 + 4*VA*(x2*x2*x2*x2)*x4 + 4*V2*(x2*x2*x2)*x3*x4 + 2*V2*(x2*x2*x2)*(x4*x4) - 4*VA*(x2*x2*x2)*(x4*x4) + 4*V2*(x1*x1)*(x2*x2)*x5 + 8*VA*(x1*x1)*(x2*x2)*x5 - 6*V2*x1*(x2*x2*x2)*x5 - 12*VA*x1*(x2*x2*x2)*x5 - 8*V2*(x2*x2*x2*x2)*x5 - 8*VA*(x2*x2*x2*x2)*x5 + 8*V2*x1*(x2*x2)*x3*x5 + 16*VA*x1*(x2*x2)*x3*x5 + 6*V2*(x2*x2*x2)*x3*x5 + 12*VA*(x2*x2*x2)*x3*x5 + 8*V2*x1*(x2*x2)*x4*x5 + 16*VA*x1*(x2*x2)*x4*x5 + 6*V2*(x2*x2*x2)*x4*x5 - 4*VA*(x2*x2*x2)*x4*x5 + 2*V2*(x2*x2)*x3*x4*x5 + 4*VA*(x2*x2)*x3*x4*x5 - 2*V2*(x2*x2)*(x4*x4)*x5 + 4*VA*(x2*x2)*(x4*x4)*x5 - 3*V2*(x1*x1)*x2*(x5*x5) - 6*VA*(x1*x1)*x2*(x5*x5) + 12*V2*x1*(x2*x2)*(x5*x5) + 24*VA*x1*(x2*x2)*(x5*x5) + 8*V2*(x2*x2*x2)*(x5*x5) + 16*VA*(x2*x2*x2)*(x5*x5) - 6*V2*x1*x2*x3*(x5*x5) - 12*VA*x1*x2*x3*(x5*x5) - 4*V2*(x2*x2)*x3*(x5*x5) - 8*VA*(x2*x2)*x3*(x5*x5) - 12*V2*x1*x2*x4*(x5*x5) - 24*VA*x1*x2*x4*(x5*x5) - 6*V2*x1*x2*(x5*x5*x5) - 12*VA*x1*x2*(x5*x5*x5) - 4*V2*(x2*x2)*(x5*x5*x5) - 8*VA*(x2*x2)*(x5*x5*x5) + 6*V2*x1*x4*(x5*x5*x5) + 12*VA*x1*x4*(x5*x5*x5) + 4*V2*(x1*x1)*(x2*x2)*x6 - 8*VA*(x1*x1)*(x2*x2)*x6 - 6*V2*x1*(x2*x2*x2)*x6 + 12*VA*x1*(x2*x2*x2)*x6 - 8*V2*(x2*x2*x2*x2)*x6 + 8*VA*(x2*x2*x2*x2)*x6 + 8*V2*x1*(x2*x2)*x3*x6 - 16*VA*x1*(x2*x2)*x3*x6 + 6*V2*(x2*x2*x2)*x3*x6 + 4*VA*(x2*x2*x2)*x3*x6 - 2*V2*(x2*x2)*(x3*x3)*x6 - 4*VA*(x2*x2)*(x3*x3)*x6 + 8*V2*x1*(x2*x2)*x4*x6 - 16*VA*x1*(x2*x2)*x4*x6 + 6*V2*(x2*x2*x2)*x4*x6 - 12*VA*(x2*x2*x2)*x4*x6 + 2*V2*(x2*x2)*x3*x4*x6 - 4*VA*(x2*x2)*x3*x4*x6 + 16*V2*x1*(x2*x2)*x5*x6 + 8*V2*(x2*x2*x2)*x5*x6 - 6*V2*x1*x2*x3*x5*x6 - 12*VA*x1*x2*x3*x5*x6 - 4*V2*(x2*x2)*x3*x5*x6 - 8*VA*(x2*x2)*x3*x5*x6 - 6*V2*x1*x2*x4*x5*x6 + 12*VA*x1*x2*x4*x5*x6 - 4*V2*(x2*x2)*x4*x5*x6 + 8*VA*(x2*x2)*x4*x5*x6 - 18*V2*x1*x2*(x5*x5)*x6 - 36*VA*x1*x2*(x5*x5)*x6 - 4*V2*(x2*x2)*(x5*x5)*x6 - 8*VA*(x2*x2)*(x5*x5)*x6 + 6*V2*x1*x3*(x5*x5)*x6 + 12*VA*x1*x3*(x5*x5)*x6 + 12*V2*x1*(x5*x5*x5)*x6 + 24*VA*x1*(x5*x5*x5)*x6 - 3*V2*(x1*x1)*x2*(x6*x6) + 6*VA*(x1*x1)*x2*(x6*x6) + 12*V2*x1*(x2*x2)*(x6*x6) - 24*VA*x1*(x2*x2)*(x6*x6) + 8*V2*(x2*x2*x2)*(x6*x6) - 16*VA*(x2*x2*x2)*(x6*x6) - 12*V2*x1*x2*x3*(x6*x6) + 24*VA*x1*x2*x3*(x6*x6) - 6*V2*x1*x2*x4*(x6*x6) + 12*VA*x1*x2*x4*(x6*x6) - 4*V2*(x2*x2)*x4*(x6*x6) + 8*VA*(x2*x2)*x4*(x6*x6) - 18*V2*x1*x2*x5*(x6*x6) + 36*VA*x1*x2*x5*(x6*x6) - 4*V2*(x2*x2)*x5*(x6*x6) + 8*VA*(x2*x2)*x5*(x6*x6) + 6*V2*x1*x4*x5*(x6*x6) - 12*VA*x1*x4*x5*(x6*x6) - 6*V2*x1*x2*(x6*x6*x6) + 12*VA*x1*x2*(x6*x6*x6) - 4*V2*(x2*x2)*(x6*x6*x6) + 8*VA*(x2*x2)*(x6*x6*x6) + 6*V2*x1*x3*(x6*x6*x6) - 12*VA*x1*x3*(x6*x6*x6) + 12*V2*x1*x5*(x6*x6*x6) - 24*VA*x1*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1)*ml2*V2*(x2*x2)*(-x2 + x5 + x6) + ml2*ml2*ml2*ml2*(x2*x2)*(2*VA*(x2 - x5 - x6) + V2*(-x5 + x6)) + ml1*ml1*ml1*ml1*(x2*x2)*(V2*(x5 - x6) + 2*VA*(-x2 + x5 + x6)) - 2*ml1*ml2*V2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6)))) - ml2*ml2*(V2*(-2*(x2*x2*x2*x2) + x2*x2*x2*(x3 + x4 + 4*x5) + x2*x2*(-2*x5*(x3 - x4 + 2*x5) + x1*(3*x5 - x6)) + 6*x1*x5*(x5 - x6)*x6 + 3*x1*x2*(-(x5*x5) + x6*x6)) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x3 + 2*x5)*(x5 + x6) + x2*(3*x3 + x4 + 6*x5 + 2*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)))) + ml1*ml1*(V2*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x5 + x6)*(x4 + 2*x6) + x2*(x3 + 3*x4 + 2*x5 + 6*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6))))) + x2*x2*(4*(Mn*Mn)*V2*x1*(x2*x2*x2*x2) - 8*(Mn*Mn)*V2*(x2*x2*x2*x2*x2) - 2*(Mn*Mn)*V2*x1*(x2*x2*x2)*x3 - 4*(Mn*Mn)*VA*x1*(x2*x2*x2)*x3 + 12*(Mn*Mn)*V2*(x2*x2*x2*x2)*x3 + 8*(Mn*Mn)*VA*(x2*x2*x2*x2)*x3 - 4*(Mn*Mn)*V2*(x2*x2*x2)*(x3*x3) - 8*(Mn*Mn)*VA*(x2*x2*x2)*(x3*x3) - 2*(Mn*Mn)*V2*x1*(x2*x2*x2)*x4 + 4*(Mn*Mn)*VA*x1*(x2*x2*x2)*x4 + 12*(Mn*Mn)*V2*(x2*x2*x2*x2)*x4 - 8*(Mn*Mn)*VA*(x2*x2*x2*x2)*x4 - 8*(Mn*Mn)*V2*(x2*x2*x2)*x3*x4 - 4*(Mn*Mn)*V2*(x2*x2*x2)*(x4*x4) + 8*(Mn*Mn)*VA*(x2*x2*x2)*(x4*x4) - 4*(Mn*Mn)*V2*(x1*x1)*(x2*x2)*x5 - 8*(Mn*Mn)*VA*(x1*x1)*(x2*x2)*x5 + V2*(x1*x1*x1)*(x2*x2)*x5 + 2*VA*(x1*x1*x1)*(x2*x2)*x5 + 4*(Mn*Mn)*V2*x1*(x2*x2*x2)*x5 + 8*(Mn*Mn)*VA*x1*(x2*x2*x2)*x5 - 2*V2*(x1*x1)*(x2*x2*x2)*x5 - 4*VA*(x1*x1)*(x2*x2*x2)*x5 + 16*(Mn*Mn)*V2*(x2*x2*x2*x2)*x5 + 16*(Mn*Mn)*VA*(x2*x2*x2*x2)*x5 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x3*x5 - 16*(Mn*Mn)*VA*x1*(x2*x2)*x3*x5 + 2*V2*(x1*x1)*(x2*x2)*x3*x5 + 4*VA*(x1*x1)*(x2*x2)*x3*x5 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x3*x5 - 24*(Mn*Mn)*VA*(x2*x2*x2)*x3*x5 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x4*x5 - 16*(Mn*Mn)*VA*x1*(x2*x2)*x4*x5 + 2*V2*(x1*x1)*(x2*x2)*x4*x5 + 4*VA*(x1*x1)*(x2*x2)*x4*x5 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x4*x5 + 8*(Mn*Mn)*VA*(x2*x2*x2)*x4*x5 - 4*(Mn*Mn)*V2*(x2*x2)*x3*x4*x5 - 8*(Mn*Mn)*VA*(x2*x2)*x3*x4*x5 + 4*(Mn*Mn)*V2*(x2*x2)*(x4*x4)*x5 - 8*(Mn*Mn)*VA*(x2*x2)*(x4*x4)*x5 + 2*(Mn*Mn)*V2*(x1*x1)*x2*(x5*x5) + 4*(Mn*Mn)*VA*(x1*x1)*x2*(x5*x5) - V2*(x1*x1*x1)*x2*(x5*x5) - 2*VA*(x1*x1*x1)*x2*(x5*x5) - 8*(Mn*Mn)*V2*x1*(x2*x2)*(x5*x5) - 16*(Mn*Mn)*VA*x1*(x2*x2)*(x5*x5) + 4*V2*(x1*x1)*(x2*x2)*(x5*x5) + 8*VA*(x1*x1)*(x2*x2)*(x5*x5) - 16*(Mn*Mn)*V2*(x2*x2*x2)*(x5*x5) - 32*(Mn*Mn)*VA*(x2*x2*x2)*(x5*x5) + 4*(Mn*Mn)*V2*x1*x2*x3*(x5*x5) + 8*(Mn*Mn)*VA*x1*x2*x3*(x5*x5) - 2*V2*(x1*x1)*x2*x3*(x5*x5) - 4*VA*(x1*x1)*x2*x3*(x5*x5) + 8*(Mn*Mn)*V2*(x2*x2)*x3*(x5*x5) + 16*(Mn*Mn)*VA*(x2*x2)*x3*(x5*x5) + 8*(Mn*Mn)*V2*x1*x2*x4*(x5*x5) + 16*(Mn*Mn)*VA*x1*x2*x4*(x5*x5) - 4*V2*(x1*x1)*x2*x4*(x5*x5) - 8*VA*(x1*x1)*x2*x4*(x5*x5) + 4*(Mn*Mn)*V2*x1*x2*(x5*x5*x5) + 8*(Mn*Mn)*VA*x1*x2*(x5*x5*x5) - 2*V2*(x1*x1)*x2*(x5*x5*x5) - 4*VA*(x1*x1)*x2*(x5*x5*x5) + 8*(Mn*Mn)*V2*(x2*x2)*(x5*x5*x5) + 16*(Mn*Mn)*VA*(x2*x2)*(x5*x5*x5) - 4*(Mn*Mn)*V2*x1*x4*(x5*x5*x5) - 8*(Mn*Mn)*VA*x1*x4*(x5*x5*x5) + 2*V2*(x1*x1)*x4*(x5*x5*x5) + 4*VA*(x1*x1)*x4*(x5*x5*x5) + 4*(ml1*ml1*ml1)*ml2*(Mn*Mn)*V2*(x2*x2)*(x2 - x5 - x6) - 4*(Mn*Mn)*V2*(x1*x1)*(x2*x2)*x6 + 8*(Mn*Mn)*VA*(x1*x1)*(x2*x2)*x6 + V2*(x1*x1*x1)*(x2*x2)*x6 - 2*VA*(x1*x1*x1)*(x2*x2)*x6 + 4*(Mn*Mn)*V2*x1*(x2*x2*x2)*x6 - 8*(Mn*Mn)*VA*x1*(x2*x2*x2)*x6 - 2*V2*(x1*x1)*(x2*x2*x2)*x6 + 4*VA*(x1*x1)*(x2*x2*x2)*x6 + 16*(Mn*Mn)*V2*(x2*x2*x2*x2)*x6 - 16*(Mn*Mn)*VA*(x2*x2*x2*x2)*x6 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x3*x6 + 16*(Mn*Mn)*VA*x1*(x2*x2)*x3*x6 + 2*V2*(x1*x1)*(x2*x2)*x3*x6 - 4*VA*(x1*x1)*(x2*x2)*x3*x6 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x3*x6 - 8*(Mn*Mn)*VA*(x2*x2*x2)*x3*x6 + 4*(Mn*Mn)*V2*(x2*x2)*(x3*x3)*x6 + 8*(Mn*Mn)*VA*(x2*x2)*(x3*x3)*x6 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x4*x6 + 16*(Mn*Mn)*VA*x1*(x2*x2)*x4*x6 + 2*V2*(x1*x1)*(x2*x2)*x4*x6 - 4*VA*(x1*x1)*(x2*x2)*x4*x6 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x4*x6 + 24*(Mn*Mn)*VA*(x2*x2*x2)*x4*x6 - 4*(Mn*Mn)*V2*(x2*x2)*x3*x4*x6 + 8*(Mn*Mn)*VA*(x2*x2)*x3*x4*x6 - 16*(Mn*Mn)*V2*x1*(x2*x2)*x5*x6 + 4*V2*(x1*x1)*(x2*x2)*x5*x6 - 16*(Mn*Mn)*V2*(x2*x2*x2)*x5*x6 + 4*(Mn*Mn)*V2*x1*x2*x3*x5*x6 + 8*(Mn*Mn)*VA*x1*x2*x3*x5*x6 - 2*V2*(x1*x1)*x2*x3*x5*x6 - 4*VA*(x1*x1)*x2*x3*x5*x6 + 8*(Mn*Mn)*V2*(x2*x2)*x3*x5*x6 + 16*(Mn*Mn)*VA*(x2*x2)*x3*x5*x6 + 4*(Mn*Mn)*V2*x1*x2*x4*x5*x6 - 8*(Mn*Mn)*VA*x1*x2*x4*x5*x6 - 2*V2*(x1*x1)*x2*x4*x5*x6 + 4*VA*(x1*x1)*x2*x4*x5*x6 + 8*(Mn*Mn)*V2*(x2*x2)*x4*x5*x6 - 16*(Mn*Mn)*VA*(x2*x2)*x4*x5*x6 + 12*(Mn*Mn)*V2*x1*x2*(x5*x5)*x6 + 24*(Mn*Mn)*VA*x1*x2*(x5*x5)*x6 - 6*V2*(x1*x1)*x2*(x5*x5)*x6 - 12*VA*(x1*x1)*x2*(x5*x5)*x6 + 8*(Mn*Mn)*V2*(x2*x2)*(x5*x5)*x6 + 16*(Mn*Mn)*VA*(x2*x2)*(x5*x5)*x6 - 4*(Mn*Mn)*V2*x1*x3*(x5*x5)*x6 - 8*(Mn*Mn)*VA*x1*x3*(x5*x5)*x6 + 2*V2*(x1*x1)*x3*(x5*x5)*x6 + 4*VA*(x1*x1)*x3*(x5*x5)*x6 - 8*(Mn*Mn)*V2*x1*(x5*x5*x5)*x6 - 16*(Mn*Mn)*VA*x1*(x5*x5*x5)*x6 + 4*V2*(x1*x1)*(x5*x5*x5)*x6 + 8*VA*(x1*x1)*(x5*x5*x5)*x6 + 2*(Mn*Mn)*V2*(x1*x1)*x2*(x6*x6) - 4*(Mn*Mn)*VA*(x1*x1)*x2*(x6*x6) - V2*(x1*x1*x1)*x2*(x6*x6) + 2*VA*(x1*x1*x1)*x2*(x6*x6) - 8*(Mn*Mn)*V2*x1*(x2*x2)*(x6*x6) + 16*(Mn*Mn)*VA*x1*(x2*x2)*(x6*x6) + 4*V2*(x1*x1)*(x2*x2)*(x6*x6) - 8*VA*(x1*x1)*(x2*x2)*(x6*x6) - 16*(Mn*Mn)*V2*(x2*x2*x2)*(x6*x6) + 32*(Mn*Mn)*VA*(x2*x2*x2)*(x6*x6) + 8*(Mn*Mn)*V2*x1*x2*x3*(x6*x6) - 16*(Mn*Mn)*VA*x1*x2*x3*(x6*x6) - 4*V2*(x1*x1)*x2*x3*(x6*x6) + 8*VA*(x1*x1)*x2*x3*(x6*x6) + 4*(Mn*Mn)*V2*x1*x2*x4*(x6*x6) - 8*(Mn*Mn)*VA*x1*x2*x4*(x6*x6) - 2*V2*(x1*x1)*x2*x4*(x6*x6) + 4*VA*(x1*x1)*x2*x4*(x6*x6) + 8*(Mn*Mn)*V2*(x2*x2)*x4*(x6*x6) - 16*(Mn*Mn)*VA*(x2*x2)*x4*(x6*x6) + 12*(Mn*Mn)*V2*x1*x2*x5*(x6*x6) - 24*(Mn*Mn)*VA*x1*x2*x5*(x6*x6) - 6*V2*(x1*x1)*x2*x5*(x6*x6) + 12*VA*(x1*x1)*x2*x5*(x6*x6) + 8*(Mn*Mn)*V2*(x2*x2)*x5*(x6*x6) - 16*(Mn*Mn)*VA*(x2*x2)*x5*(x6*x6) - 4*(Mn*Mn)*V2*x1*x4*x5*(x6*x6) + 8*(Mn*Mn)*VA*x1*x4*x5*(x6*x6) + 2*V2*(x1*x1)*x4*x5*(x6*x6) - 4*VA*(x1*x1)*x4*x5*(x6*x6) + 4*(Mn*Mn)*V2*x1*x2*(x6*x6*x6) - 8*(Mn*Mn)*VA*x1*x2*(x6*x6*x6) - 2*V2*(x1*x1)*x2*(x6*x6*x6) + 4*VA*(x1*x1)*x2*(x6*x6*x6) + 8*(Mn*Mn)*V2*(x2*x2)*(x6*x6*x6) - 16*(Mn*Mn)*VA*(x2*x2)*(x6*x6*x6) - 4*(Mn*Mn)*V2*x1*x3*(x6*x6*x6) + 8*(Mn*Mn)*VA*x1*x3*(x6*x6*x6) + 2*V2*(x1*x1)*x3*(x6*x6*x6) - 4*VA*(x1*x1)*x3*(x6*x6*x6) - 8*(Mn*Mn)*V2*x1*x5*(x6*x6*x6) + 16*(Mn*Mn)*VA*x1*x5*(x6*x6*x6) + 4*V2*(x1*x1)*x5*(x6*x6*x6) - 8*VA*(x1*x1)*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1*ml1)*(Mn*Mn)*(x2*x2)*(2*VA*(x2 - x5 - x6) + V2*(-x5 + x6)) + 2*(ml2*ml2*ml2*ml2)*(Mn*Mn)*(x2*x2)*(V2*(x5 - x6) + 2*VA*(-x2 + x5 + x6)) + 2*ml1*ml2*V2*(2*(ml2*ml2)*(Mn*Mn)*(x2*x2)*(x2 - x5 - x6) - x1*x1*(x2*x2*x2 - 2*(x2*x2)*(x5 + x6) - 2*x5*x6*(x5 + x6) + x2*(x5*x5 + 4*x5*x6 + x6*x6)) + 2*(Mn*Mn)*(-(x1*(x2*x2*x2 + 2*(x2*x2)*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6))))) + ml2*ml2*(x1*x1*(-(V2*(x5 - x6)*(x2*x2 + 2*x5*x6 - x2*(x5 + x6))) - 2*VA*(x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) + 2*(Mn*Mn)*(V2*(-2*(x2*x2*x2*x2) + x2*x2*x2*(x3 + x4 + 4*x5) + 2*x1*x5*(x5 - x6)*x6 + x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(-2*x5*(x3 - x4 + 2*x5) + x1*(x5 + x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x3 + 2*x5)*(x5 + x6) + x2*(3*x3 + x4 + 6*x5 + 2*x6)) + x1*(x2*x2*x2 + x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))))) + ml1*ml1*(x1*x1*(V2*(x5 - x6)*(x2*x2 + 2*x5*x6 - x2*(x5 + x6)) + 2*VA*(x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) - 2*(Mn*Mn)*(V2*(2*(x2*x2*x2*x2) + 2*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + x1*x2*(-(x5*x5) + x6*x6) - x2*x2*(2*(x3 - x4 - 2*x6)*x6 + x1*(x5 + x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x5 + x6)*(x4 + 2*x6) + x2*(x3 + 3*x4 + 2*x5 + 6*x6)) + x1*(x2*x2*x2 + x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6)))))) + A2*(2*(Enu*Enu)*(Mn*Mn)*x1*(-2*x1*(x2*x2*x2*x2) + 4*(x2*x2*x2*x2*x2) + x1*(x2*x2*x2)*x3 - 6*(x2*x2*x2*x2)*x3 + 2*(x2*x2*x2)*(x3*x3) + x1*(x2*x2*x2)*x4 - 6*(x2*x2*x2*x2)*x4 + 4*(x2*x2*x2)*x3*x4 + 2*(x2*x2*x2)*(x4*x4) + 4*(x1*x1)*(x2*x2)*x5 - 6*x1*(x2*x2*x2)*x5 - 8*(x2*x2*x2*x2)*x5 + 8*x1*(x2*x2)*x3*x5 + 6*(x2*x2*x2)*x3*x5 + 8*x1*(x2*x2)*x4*x5 + 6*(x2*x2*x2)*x4*x5 + 2*(x2*x2)*x3*x4*x5 - 2*(x2*x2)*(x4*x4)*x5 - 3*(x1*x1)*x2*(x5*x5) + 12*x1*(x2*x2)*(x5*x5) + 8*(x2*x2*x2)*(x5*x5) - 6*x1*x2*x3*(x5*x5) - 4*(x2*x2)*x3*(x5*x5) - 12*x1*x2*x4*(x5*x5) - 6*x1*x2*(x5*x5*x5) - 4*(x2*x2)*(x5*x5*x5) + 6*x1*x4*(x5*x5*x5) + 2*(ml1*ml1*ml1)*ml2*(x2*x2)*(x2 - x5 - x6) + ml1*ml1*ml1*ml1*(x2*x2)*(x5 - x6) + 4*(x1*x1)*(x2*x2)*x6 - 6*x1*(x2*x2*x2)*x6 - 8*(x2*x2*x2*x2)*x6 + 8*x1*(x2*x2)*x3*x6 + 6*(x2*x2*x2)*x3*x6 - 2*(x2*x2)*(x3*x3)*x6 + 8*x1*(x2*x2)*x4*x6 + 6*(x2*x2*x2)*x4*x6 + 2*(x2*x2)*x3*x4*x6 + 16*x1*(x2*x2)*x5*x6 + 8*(x2*x2*x2)*x5*x6 - 6*x1*x2*x3*x5*x6 - 4*(x2*x2)*x3*x5*x6 - 6*x1*x2*x4*x5*x6 - 4*(x2*x2)*x4*x5*x6 - 18*x1*x2*(x5*x5)*x6 - 4*(x2*x2)*(x5*x5)*x6 + 6*x1*x3*(x5*x5)*x6 + 12*x1*(x5*x5*x5)*x6 - 3*(x1*x1)*x2*(x6*x6) + 12*x1*(x2*x2)*(x6*x6) + 8*(x2*x2*x2)*(x6*x6) - 12*x1*x2*x3*(x6*x6) - 6*x1*x2*x4*(x6*x6) - 4*(x2*x2)*x4*(x6*x6) - 18*x1*x2*x5*(x6*x6) - 4*(x2*x2)*x5*(x6*x6) + 6*x1*x4*x5*(x6*x6) - 6*x1*x2*(x6*x6*x6) - 4*(x2*x2)*(x6*x6*x6) + 6*x1*x3*(x6*x6*x6) + 12*x1*x5*(x6*x6*x6) + ml2*ml2*ml2*ml2*(x2*x2)*(-x5 + x6) + ml2*ml2*(2*(x2*x2*x2*x2) - x2*x2*x2*(x3 + x4 + 4*x5) + 6*x1*x5*x6*(-x5 + x6) + 3*x1*x2*(x5*x5 - x6*x6) + x2*x2*(2*x5*(x3 - x4 + 2*x5) + x1*(-3*x5 + x6))) + ml1*ml1*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*ml1*ml2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6))))) - 2*Enu*Mn*x1*x2*(-2*x1*(x2*x2*x2*x2) + 4*(x2*x2*x2*x2*x2) + x1*(x2*x2*x2)*x3 - 6*(x2*x2*x2*x2)*x3 + 2*(x2*x2*x2)*(x3*x3) + x1*(x2*x2*x2)*x4 - 6*(x2*x2*x2*x2)*x4 + 4*(x2*x2*x2)*x3*x4 + 2*(x2*x2*x2)*(x4*x4) + 4*(x1*x1)*(x2*x2)*x5 - 6*x1*(x2*x2*x2)*x5 - 8*(x2*x2*x2*x2)*x5 + 8*x1*(x2*x2)*x3*x5 + 6*(x2*x2*x2)*x3*x5 + 8*x1*(x2*x2)*x4*x5 + 6*(x2*x2*x2)*x4*x5 + 2*(x2*x2)*x3*x4*x5 - 2*(x2*x2)*(x4*x4)*x5 - 3*(x1*x1)*x2*(x5*x5) + 12*x1*(x2*x2)*(x5*x5) + 8*(x2*x2*x2)*(x5*x5) - 6*x1*x2*x3*(x5*x5) - 4*(x2*x2)*x3*(x5*x5) - 12*x1*x2*x4*(x5*x5) - 6*x1*x2*(x5*x5*x5) - 4*(x2*x2)*(x5*x5*x5) + 6*x1*x4*(x5*x5*x5) + 2*(ml1*ml1*ml1)*ml2*(x2*x2)*(x2 - x5 - x6) + ml1*ml1*ml1*ml1*(x2*x2)*(x5 - x6) + 4*(x1*x1)*(x2*x2)*x6 - 6*x1*(x2*x2*x2)*x6 - 8*(x2*x2*x2*x2)*x6 + 8*x1*(x2*x2)*x3*x6 + 6*(x2*x2*x2)*x3*x6 - 2*(x2*x2)*(x3*x3)*x6 + 8*x1*(x2*x2)*x4*x6 + 6*(x2*x2*x2)*x4*x6 + 2*(x2*x2)*x3*x4*x6 + 16*x1*(x2*x2)*x5*x6 + 8*(x2*x2*x2)*x5*x6 - 6*x1*x2*x3*x5*x6 - 4*(x2*x2)*x3*x5*x6 - 6*x1*x2*x4*x5*x6 - 4*(x2*x2)*x4*x5*x6 - 18*x1*x2*(x5*x5)*x6 - 4*(x2*x2)*(x5*x5)*x6 + 6*x1*x3*(x5*x5)*x6 + 12*x1*(x5*x5*x5)*x6 - 3*(x1*x1)*x2*(x6*x6) + 12*x1*(x2*x2)*(x6*x6) + 8*(x2*x2*x2)*(x6*x6) - 12*x1*x2*x3*(x6*x6) - 6*x1*x2*x4*(x6*x6) - 4*(x2*x2)*x4*(x6*x6) - 18*x1*x2*x5*(x6*x6) - 4*(x2*x2)*x5*(x6*x6) + 6*x1*x4*x5*(x6*x6) - 6*x1*x2*(x6*x6*x6) - 4*(x2*x2)*(x6*x6*x6) + 6*x1*x3*(x6*x6*x6) + 12*x1*x5*(x6*x6*x6) + ml2*ml2*ml2*ml2*(x2*x2)*(-x5 + x6) + ml2*ml2*(2*(x2*x2*x2*x2) - x2*x2*x2*(x3 + x4 + 4*x5) + 6*x1*x5*x6*(-x5 + x6) + 3*x1*x2*(x5*x5 - x6*x6) + x2*x2*(2*x5*(x3 - x4 + 2*x5) + x1*(-3*x5 + x6))) + ml1*ml1*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*ml1*ml2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6))))) + x2*x2*(-4*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2*x2) + 4*(Mn*Mn)*x1*(x2*x2*x2*x2) - 8*(Mn*Mn)*(x2*x2*x2*x2*x2) + 2*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2)*x3 - 2*(Mn*Mn)*x1*(x2*x2*x2)*x3 + 12*(Mn*Mn)*(x2*x2*x2*x2)*x3 - 4*(Mn*Mn)*(x2*x2*x2)*(x3*x3) + 2*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2)*x4 - 2*(Mn*Mn)*x1*(x2*x2*x2)*x4 + 12*(Mn*Mn)*(x2*x2*x2*x2)*x4 - 8*(Mn*Mn)*(x2*x2*x2)*x3*x4 - 4*(Mn*Mn)*(x2*x2*x2)*(x4*x4) + 2*(ml2*ml2*ml2*ml2)*(Mn*Mn)*(x2*x2)*x5 + 2*(ml2*ml2)*(Mn*Mn)*x1*(x2*x2)*x5 - ml2*ml2*(x1*x1)*(x2*x2)*x5 - 4*(Mn*Mn)*(x1*x1)*(x2*x2)*x5 + x1*x1*x1*(x2*x2)*x5 + 8*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2)*x5 + 4*(Mn*Mn)*x1*(x2*x2*x2)*x5 - 2*(x1*x1)*(x2*x2*x2)*x5 + 16*(Mn*Mn)*(x2*x2*x2*x2)*x5 - 4*(ml2*ml2)*(Mn*Mn)*(x2*x2)*x3*x5 - 8*(Mn*Mn)*x1*(x2*x2)*x3*x5 + 2*(x1*x1)*(x2*x2)*x3*x5 - 12*(Mn*Mn)*(x2*x2*x2)*x3*x5 + 4*(ml2*ml2)*(Mn*Mn)*(x2*x2)*x4*x5 - 8*(Mn*Mn)*x1*(x2*x2)*x4*x5 + 2*(x1*x1)*(x2*x2)*x4*x5 - 12*(Mn*Mn)*(x2*x2*x2)*x4*x5 - 4*(Mn*Mn)*(x2*x2)*x3*x4*x5 + 4*(Mn*Mn)*(x2*x2)*(x4*x4)*x5 - 2*(ml2*ml2)*(Mn*Mn)*x1*x2*(x5*x5) + ml2*ml2*(x1*x1)*x2*(x5*x5) + 2*(Mn*Mn)*(x1*x1)*x2*(x5*x5) - x1*x1*x1*x2*(x5*x5) - 8*(ml2*ml2)*(Mn*Mn)*(x2*x2)*(x5*x5) - 8*(Mn*Mn)*x1*(x2*x2)*(x5*x5) + 4*(x1*x1)*(x2*x2)*(x5*x5) - 16*(Mn*Mn)*(x2*x2*x2)*(x5*x5) + 4*(Mn*Mn)*x1*x2*x3*(x5*x5) - 2*(x1*x1)*x2*x3*(x5*x5) + 8*(Mn*Mn)*(x2*x2)*x3*(x5*x5) + 8*(Mn*Mn)*x1*x2*x4*(x5*x5) - 4*(x1*x1)*x2*x4*(x5*x5) + 4*(Mn*Mn)*x1*x2*(x5*x5*x5) - 2*(x1*x1)*x2*(x5*x5*x5) + 8*(Mn*Mn)*(x2*x2)*(x5*x5*x5) - 4*(Mn*Mn)*x1*x4*(x5*x5*x5) + 2*(x1*x1)*x4*(x5*x5*x5) - 2*(ml2*ml2*ml2*ml2)*(Mn*Mn)*(x2*x2)*x6 + 2*(ml2*ml2)*(Mn*Mn)*x1*(x2*x2)*x6 + ml2*ml2*(x1*x1)*(x2*x2)*x6 - 4*(Mn*Mn)*(x1*x1)*(x2*x2)*x6 + x1*x1*x1*(x2*x2)*x6 + 4*(Mn*Mn)*x1*(x2*x2*x2)*x6 - 2*(x1*x1)*(x2*x2*x2)*x6 + 16*(Mn*Mn)*(x2*x2*x2*x2)*x6 - 8*(Mn*Mn)*x1*(x2*x2)*x3*x6 + 2*(x1*x1)*(x2*x2)*x3*x6 - 12*(Mn*Mn)*(x2*x2*x2)*x3*x6 + 4*(Mn*Mn)*(x2*x2)*(x3*x3)*x6 - 8*(Mn*Mn)*x1*(x2*x2)*x4*x6 + 2*(x1*x1)*(x2*x2)*x4*x6 - 12*(Mn*Mn)*(x2*x2*x2)*x4*x6 - 4*(Mn*Mn)*(x2*x2)*x3*x4*x6 - 16*(Mn*Mn)*x1*(x2*x2)*x5*x6 + 4*(x1*x1)*(x2*x2)*x5*x6 - 16*(Mn*Mn)*(x2*x2*x2)*x5*x6 + 4*(Mn*Mn)*x1*x2*x3*x5*x6 - 2*(x1*x1)*x2*x3*x5*x6 + 8*(Mn*Mn)*(x2*x2)*x3*x5*x6 + 4*(Mn*Mn)*x1*x2*x4*x5*x6 - 2*(x1*x1)*x2*x4*x5*x6 + 8*(Mn*Mn)*(x2*x2)*x4*x5*x6 + 4*(ml2*ml2)*(Mn*Mn)*x1*(x5*x5)*x6 - 2*(ml2*ml2)*(x1*x1)*(x5*x5)*x6 + 12*(Mn*Mn)*x1*x2*(x5*x5)*x6 - 6*(x1*x1)*x2*(x5*x5)*x6 + 8*(Mn*Mn)*(x2*x2)*(x5*x5)*x6 - 4*(Mn*Mn)*x1*x3*(x5*x5)*x6 + 2*(x1*x1)*x3*(x5*x5)*x6 - 8*(Mn*Mn)*x1*(x5*x5*x5)*x6 + 4*(x1*x1)*(x5*x5*x5)*x6 + 2*(ml2*ml2)*(Mn*Mn)*x1*x2*(x6*x6) - ml2*ml2*(x1*x1)*x2*(x6*x6) + 2*(Mn*Mn)*(x1*x1)*x2*(x6*x6) - x1*x1*x1*x2*(x6*x6) - 8*(Mn*Mn)*x1*(x2*x2)*(x6*x6) + 4*(x1*x1)*(x2*x2)*(x6*x6) - 16*(Mn*Mn)*(x2*x2*x2)*(x6*x6) + 8*(Mn*Mn)*x1*x2*x3*(x6*x6) - 4*(x1*x1)*x2*x3*(x6*x6) + 4*(Mn*Mn)*x1*x2*x4*(x6*x6) - 2*(x1*x1)*x2*x4*(x6*x6) + 8*(Mn*Mn)*(x2*x2)*x4*(x6*x6) - 4*(ml2*ml2)*(Mn*Mn)*x1*x5*(x6*x6) + 2*(ml2*ml2)*(x1*x1)*x5*(x6*x6) + 12*(Mn*Mn)*x1*x2*x5*(x6*x6) - 6*(x1*x1)*x2*x5*(x6*x6) + 8*(Mn*Mn)*(x2*x2)*x5*(x6*x6) - 4*(Mn*Mn)*x1*x4*x5*(x6*x6) + 2*(x1*x1)*x4*x5*(x6*x6) + 4*(Mn*Mn)*x1*x2*(x6*x6*x6) - 2*(x1*x1)*x2*(x6*x6*x6) + 8*(Mn*Mn)*(x2*x2)*(x6*x6*x6) - 4*(Mn*Mn)*x1*x3*(x6*x6*x6) + 2*(x1*x1)*x3*(x6*x6*x6) - 8*(Mn*Mn)*x1*x5*(x6*x6*x6) + 4*(x1*x1)*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1*ml1)*(Mn*Mn)*(x2*x2)*(-x5 + x6) + 4*(ml1*ml1*ml1)*ml2*(Mn*Mn)*(x2*x2)*(-x2 + x5 + x6) + ml1*ml1*(x1*x1*(x5 - x6)*(x2*x2 + 2*x5*x6 - x2*(x5 + x6)) + 2*(Mn*Mn)*(-2*(x2*x2*x2*x2) + 2*x1*x5*x6*(-x5 + x6) + x2*x2*x2*(x3 + x4 + 4*x6) + x1*x2*(x5*x5 - x6*x6) + x2*x2*(2*(x3 - x4 - 2*x6)*x6 + x1*(x5 + x6)))) - 2*ml1*ml2*(2*(ml2*ml2)*(Mn*Mn)*(x2*x2)*(x2 - x5 - x6) - x1*x1*(x2*x2*x2 - 2*(x2*x2)*(x5 + x6) - 2*x5*x6*(x5 + x6) + x2*(x5*x5 + 4*x5*x6 + x6*x6)) + 2*(Mn*Mn)*(-(x1*(x2*x2*x2 + 2*(x2*x2)*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6)))))))))*(Z*Z))/(Enu*Enu*(Mn*Mn)*(x1*x1)*(x2*x2*x2*x2)*((x1 + 2*x3)*(x1 + 2*x3))*((x1 + 2*x4)*(x1 + 2*x4)));

  ///////////////////////////////////////////////////////////////////////////////////
  // JACOBIANS AND PREFACTORS 
  dsigma *= Jacob;

  // Changing units to zeptobarn = zb = 1e-45 cm^2
  dsigma *= (GeV2_to_cm2*1e45);

  myMC->xphys[0] = x1;
  myMC->xphys[1] = x2;
  myMC->xphys[2] = x3;
  myMC->xphys[3] = x4;
  myMC->xphys[4] = x5;
  myMC->xphys[5] = x6;
  myMC->xphys[6] = x7;
  myMC->xphys[7] = x8;
  myMC->xphys[8] = Enu;
  myMC->xphys[9] = dsigma;

  return 0;
}

int DifFromIntToPhysical_Efixed(const cubareal xx[], void *MC){

  long double u1_s    = xx[0];
  long double u2_s    = xx[1];
  long double u3_s    = xx[2];
  long double PHI2_s  = xx[3];
  long double x5_s    = xx[4];
  long double u6_s    = xx[5];
  long double x7_s    = xx[6];
  long double x8_s    = xx[7];

  // Random angles needed to specify 4 vectors
  long double x7 = x7_s*2*M_PI;
  long double x8 = x8_s*2*M_PI;

  //////////////////////////////////////////////////
  // get data from user

  tridentMC* myMC = (tridentMC *)MC;
  
  long double ml1 = myMC->ml1;
  long double ml2 = myMC->ml2;
  long double A = myMC->A;
  long double Z = myMC->Z;
  long double Mn = myMC->Mn;
  long double mzprime = myMC->mzprime;

  long double Diag11 = myMC->terms[0];
  long double Diag22 = myMC->terms[1];
  long double Diag12 = myMC->terms[2];
  long double V2     = myMC->terms[3];
  long double A2     = myMC->terms[4];
  long double VA     = myMC->terms[5];
  long double BSM    = myMC->terms[6];

  ////////////////////////////////////////////////////
  // FIXED neutrino energy 
  long double Enu = myMC->nu_energy;


  //////////////////////////////////
  long double x1_u = 2*SQR(Enu) / ( 1 + 2 * Enu/Mn ) * ( 1 - SQR(ml1+ml2)/2.0/SQR(Enu) * (1 + Enu/Mn) +
                  sqrt( SQR(1 - SQR(ml1 + ml2)/2.0/SQR(Enu)*(1+Enu/Mn) ) 
                    - pow(ml1 + ml2, 4)/4.0/pow(Enu,4)*(1 + 2*Enu/Mn) )  );
  long double x1_l = pow(ml1+ml2,4)/(x1_u)/(1.0+2*Enu/Mn);

  long double u1_l = log(x1_l); 
  long double u1_u = log(x1_u); 
  
  long double u1 = u1_s*(u1_u-u1_l) + u1_l;

  long double x1 = exp(u1);

  //////////////////////////////////
  long double x2_l = 1.0/2.0 * ( x1 + SQR(ml1 + ml2) );
  long double x2_u = Enu * (sqrt( x1 + SQR(x1)/4.0/SQR(Mn) ) - x1/2.0/Mn );
  // long double x2 = x2_s*(x2_u-x2_l) + x2_l;

  long double u2_u = (1.0+2.0*Enu/Mn)*(x1_u-x1)*(x1-x1_l)  / ( SQR(ml1 + ml2) + x1*(1+Enu/Mn) + 2*Enu*sqrt(x1+SQR(x1)/4.0/SQR(Mn)) );
  long double u2_l = 0.0;

  long double u2 = u2_s*(u2_u-u2_l) + u2_l;

  long double x2 = 0.5*(u2 + SQR(ml1 + ml2) + x1);

  /////////////////
  long double x5_l = x2/2.0/(2*x2-x1) * (  2*x2 - x1 + SQR(ml2) - SQR(ml1)  
        - sqrt( SQR( 2*x2-x1 + SQR(ml2) - SQR(ml1) ) - 4*SQR(ml2)*(2*x2-x1) ) );
  
  long double x5_u = x2/2.0/(2*x2-x1) * (  2*x2 - x1 + SQR(ml2) - SQR(ml1)  
        + sqrt( SQR( 2*x2-x1 + SQR(ml2) - SQR(ml1) ) - 4*SQR(ml2)*(2*x2-x1) ) );
  

  long double x5 = x5_s*(x5_u-x5_l) + x5_l; 

  ////////////////
  long double x3_l = SQR(ml2)/2.0 * x2/x5 - x1/2.0 * x5/x2;

  long double x3_u = 0.5*(2*x2 - x1 + SQR(ml2) - SQR(ml1) )  - x5;
  
  long double u3_l = log(2*x3_l + x1);
  long double u3_u = log(2*x3_u + x1);

  long double u3 = u3_s*(u3_u-u3_l) + u3_l;

  long double x3 = 0.5*(exp(u3)- x1);

  ////////////////


  ////////////////
  // Useful Definitions in Frame p1vec + qvec - p3vec = 0
  ////////////////
  long double Wc2 = 2*(x2-x3-x5) - x1 + SQR(ml2);

  long double E1 =  (x2 - x5)/sqrt(Wc2) ;
  long double E4 =  (Wc2+SQR(ml1))/2.0/sqrt(Wc2) ;
  long double q0 =  (x2-x1-x3)/sqrt(Wc2) ;
  long double qvec = sqrt(SQR(q0) + x1)  ;
  long double p4vec = (Wc2 - SQR(ml1))/2.0/sqrt(Wc2)  ;

  long double Cq = (q0*E1 - x2)/qvec/E1;
  long double Sq = sqrt(1 - SQR(Cq));

  long double E2 = (Wc2 - SQR(ml1))/2/sqrt(Wc2);


  ////////////////

  long double m6, u6, u6_l, u6_u, prop, jacob_u6;
  long double m6_l = 0;
  long double m6_u = 2*E1*E2;

  switch ((int) BSM) 
  {
    case (SMonly):
      u6_l = m6_l;
      u6_u = m6_u;
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6 = u6;
      jacob_u6 = 1.0;

      V2     = myMC->terms[3];
      A2     = myMC->terms[4];
      VA     = myMC->terms[5];
      break;

    case (INTERFERENCE):

      //////////////////////////////////////////////////////////
      // BSM PROPAGATOR 
      // long double prop = (1.0/(-2*(m6/mzprime/mzprime) - 1.0))/mzprime/mzprime;
      u6_l = log(1.0/(2.0*m6_u + mzprime*mzprime));
      u6_u = log(1.0/(2.0*m6_l + mzprime*mzprime));
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6   = -mzprime*mzprime/2.0 + exp(-u6)/2.0;

      prop = -exp(u6);
      jacob_u6 = exp(-u6)/2.0;

      V2     = myMC->terms[3] * prop;
      A2     = myMC->terms[4] * prop;
      VA     = myMC->terms[5] * prop;
      break;

    case (BSMonly):

      u6_l = 1.0/(2.0*m6_u + mzprime*mzprime);
      u6_u = 1.0/(2.0*m6_l + mzprime*mzprime);
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6   = -mzprime*mzprime/2.0 + 1.0/2.0/u6;

      prop = -u6;
      jacob_u6 = 1.0/(u6*u6)/2.0;

      V2     = myMC->terms[3] * prop*prop;
      A2     = myMC->terms[4] * prop*prop;
      VA     = myMC->terms[5] * prop*prop;
      break;
    
    case (SMandBSM):
     
      u6_l = m6_l;
      u6_u = m6_u;
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6 = u6;
      jacob_u6 = 1.0;

      prop = 1.0/(2.0*m6 + mzprime*mzprime);

      V2     = (myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop) ;
      A2     = (myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop) ;
      VA     = (myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop);
      break;

    default:
      std::cout << "Could not determine what contributions (BSM) to include." << std::endl;
      break;
  }
  
  long double C2 = 1.0 - m6 / E1 / E2;
  long double S2 = sqrt(1.0 - SQR(C2));


  ////////////////
  long double PHI2 = PHI2_s*(2*M_PI-0.0) + 0.0;

  long double x4 = q0*E4 + (Sq*S2*cos(PHI2) + Cq*C2)*E2*qvec;

  long double C4 = (E4*q0 - x4)/(p4vec*qvec);
  long double S4 = sqrt(1.0 - SQR(C4));

  long double x6 = E1*E4 + E1*E2*C2; 

  //////////////////////////////////////////////////////////
  // Multiply by the appropriate Jacobians for our invariants and phase space factors
  long double Jacob =  (Wc2 - SQR(ml1))/8.0/Wc2/E1/E2/x2/pow(2*M_PI,6.0);

  // Jacobian due to change from CSW(+Matheus) to better integration variables AND VEGAS
  Jacob *=  (u1_u-u1_l)*exp(u1)*
            (u2_u-u2_l)*0.5*
            (u3_u-u3_l)*0.5*exp(u3)*
            (2*M_PI)*
            (x5_u-x5_l)*
            (u6_u-u6_l)*jacob_u6*
            (2*M_PI)*
            (2*M_PI);


  //////////////////////////////////////////////////////////
  // MATRIX ELEMENT ITSELF
  long double HH1, HH2;
  if (Z == 0)
  {
    HH1 = H1_n(sqrt(x1));
    HH2 = H2_n(sqrt(x1));
  }
  else if (Z == 1)
  {
    HH1 = H1_p(sqrt(x1));
    HH2 = H2_p(sqrt(x1));
  }
  else
  {
    std::cout<<"Error! Invalid Z. Diffractive regime has Z=0 or Z=1. "<<std::endl;
  }
  // includes Z*Z for the coherent enhancement
  long double dsigma = (4*(alphaQED*alphaQED)*(Gf*Gf)*(Diag22*((x1 + 2*x4)*(x1 + 2*x4))*(2*VA*(-2*(ml2*ml2)*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) - 2*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) + 2*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5)*x5 - 4*(x1*x1)*(x2*x2)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5 + 8*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*(x5*x5) - 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*(x5*x5*x5) - 3*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 4*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x6 + 4*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(-x1 + x2 - x3 - x4)*x6 + 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x5*x5)*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 - 2*(x1*x1)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6)*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x3*x5*(x2 - x5 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x3*x5*(x2 - x5 - x6)*x6 + 2*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) - A2*(2*(ml2*ml2)*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) + 2*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) - 2*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5)*x5 + 4*(x1*x1)*(x2*x2)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5 - 8*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*(x5*x5) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*(x5*x5*x5) + 3*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 4*ml1*(ml2*ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(x2 - x5 - x6) - 4*ml1*ml2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6) + 4*ml1*ml2*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6) + 4*ml1*ml2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*(x2 - x5 - x6) - 4*ml1*ml2*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x5*x5)*(x2 - x5 - x6) + 4*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x6 + 4*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(-x1 + x2 - x3 - x4)*x6 + 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x5*x5)*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 - 2*(x1*x1)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6)*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x3*x5*(x2 - x5 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x3*x5*(x2 - x5 - x6)*x6 + 2*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) - V2*(2*(ml2*ml2)*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) + 2*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) - 2*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5)*x5 + 4*(x1*x1)*(x2*x2)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5 - 8*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*(x5*x5) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*(x5*x5*x5) + 3*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 4*ml1*(ml2*ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(x2 - x5 - x6) + 4*ml1*ml2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6) - 4*ml1*ml2*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6) - 4*ml1*ml2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*(x2 - x5 - x6) + 4*ml1*ml2*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x5*x5)*(x2 - x5 - x6) + 4*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x6 + 4*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(-x1 + x2 - x3 - x4)*x6 + 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x5*x5)*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 - 2*(x1*x1)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6)*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x3*x5*(x2 - x5 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x3*x5*(x2 - x5 - x6)*x6 + 2*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6))) + Diag11*((x1 + 2*x3)*(x1 + 2*x3))*(2*VA*(-4*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x5 - 4*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x4*x5 - 2*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*(ml1*ml1)*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) + 2*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x4*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) - 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 + 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*x6 + 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5*x6 + 2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 - 2*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6)*x6 + 4*(x1*x1)*(x2*x2)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 + 2*(x1*x1)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6)*x6 - 2*(x1*x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6)*x6 + 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x4*x5*(x2 - x5 - x6)*x6 - 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x4*x5*(x2 - x5 - x6)*x6 + 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) - 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*(x6*x6) - 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5*(x6*x6) - 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) - 8*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*(x6*x6) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*(x6*x6*x6) + 3*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) - A2*(4*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x5 + 4*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x4*x5 + 2*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*(ml1*ml1)*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) + 2*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x4*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) + 4*(ml1*ml1*ml1)*ml2*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x4*(x2 - x5 - x6) + 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5*x6 - 2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 - 2*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6)*x6 + 4*(x1*x1)*(x2*x2)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 - 4*ml1*ml2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6)*x6 + 4*ml1*ml2*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x2 - x5 - x6)*x6 - 2*(x1*x1)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6)*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x4*x5*(x2 - x5 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x4*x5*(x2 - x5 - x6)*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*(x6*x6) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5*(x6*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) - 8*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*(x6*x6) + 4*ml1*ml2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6)*(x6*x6) - 4*ml1*ml2*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x2 - x5 - x6)*(x6*x6) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*(x6*x6*x6) + 3*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) - V2*(4*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x5 + 4*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x4*x5 + 2*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*(ml1*ml1)*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) + 2*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x4*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) - 4*(ml1*ml1*ml1)*ml2*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x4*(x2 - x5 - x6) + 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5*x6 - 2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 - 2*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6)*x6 + 4*(x1*x1)*(x2*x2)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 + 4*ml1*ml2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6)*x6 - 4*ml1*ml2*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x2 - x5 - x6)*x6 - 2*(x1*x1)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6)*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x4*x5*(x2 - x5 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x4*x5*(x2 - x5 - x6)*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*(x6*x6) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5*(x6*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) - 8*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*(x6*x6) - 4*ml1*ml2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6)*(x6*x6) + 4*ml1*ml2*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x2 - x5 - x6)*(x6*x6) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*(x6*x6*x6) + 3*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6))) + 2*Diag12*(x1 + 2*x3)*(x1 + 2*x4)*(HH2*x1*(x2*x2*x2*x2)*(-2*V2*x1*(x2*x2) + 4*V2*(x2*x2*x2) + V2*x1*x2*x3 + 2*VA*x1*x2*x3 - 6*V2*(x2*x2)*x3 - 4*VA*(x2*x2)*x3 + 2*V2*x2*(x3*x3) + 4*VA*x2*(x3*x3) + V2*x1*x2*x4 - 2*VA*x1*x2*x4 - 6*V2*(x2*x2)*x4 + 4*VA*(x2*x2)*x4 + 4*V2*x2*x3*x4 + 2*V2*x2*(x4*x4) - 4*VA*x2*(x4*x4) + V2*(x1*x1)*x5 + 2*VA*(x1*x1)*x5 - 8*V2*(x2*x2)*x5 - 8*VA*(x2*x2)*x5 + 2*V2*x1*x3*x5 + 4*VA*x1*x3*x5 + 6*V2*x2*x3*x5 + 12*VA*x2*x3*x5 + 2*V2*x1*x4*x5 + 4*VA*x1*x4*x5 + 6*V2*x2*x4*x5 - 4*VA*x2*x4*x5 + 2*V2*x3*x4*x5 + 4*VA*x3*x4*x5 - 2*V2*(x4*x4)*x5 + 4*VA*(x4*x4)*x5 + 8*V2*x2*(x5*x5) + 16*VA*x2*(x5*x5) - 4*V2*x3*(x5*x5) - 8*VA*x3*(x5*x5) - 4*V2*(x5*x5*x5) - 8*VA*(x5*x5*x5) + V2*(x1*x1)*x6 - 2*VA*(x1*x1)*x6 - 8*V2*(x2*x2)*x6 + 8*VA*(x2*x2)*x6 + 2*V2*x1*x3*x6 - 4*VA*x1*x3*x6 + 6*V2*x2*x3*x6 + 4*VA*x2*x3*x6 - 2*V2*(x3*x3)*x6 - 4*VA*(x3*x3)*x6 + 2*V2*x1*x4*x6 - 4*VA*x1*x4*x6 + 6*V2*x2*x4*x6 - 12*VA*x2*x4*x6 + 2*V2*x3*x4*x6 - 4*VA*x3*x4*x6 + 4*V2*x1*x5*x6 + 8*V2*x2*x5*x6 - 4*V2*x3*x5*x6 - 8*VA*x3*x5*x6 - 4*V2*x4*x5*x6 + 8*VA*x4*x5*x6 - 4*V2*(x5*x5)*x6 - 8*VA*(x5*x5)*x6 + 8*V2*x2*(x6*x6) - 16*VA*x2*(x6*x6) - 4*V2*x4*(x6*x6) + 8*VA*x4*(x6*x6) - 4*V2*x5*(x6*x6) + 8*VA*x5*(x6*x6) - 4*V2*(x6*x6*x6) + 8*VA*(x6*x6*x6) + 2*(ml1*ml1*ml1)*ml2*V2*(-x2 + x5 + x6) - 2*ml1*ml2*V2*(-2*x1*x2 + 4*(x2*x2) - 3*x2*x3 - 3*x2*x4 - 4*x2*x5 + x3*x5 + x4*x5 + 2*(x5*x5) + ml2*ml2*(x2 - x5 - x6) - 4*x2*x6 + x3*x6 + x4*x6 + 4*x5*x6 + 2*(x6*x6)) + ml2*ml2*ml2*ml2*(2*VA*(x2 - x5 - x6) + V2*(-x5 + x6)) + ml1*ml1*ml1*ml1*(V2*(x5 - x6) + 2*VA*(-x2 + x5 + x6)) + ml1*ml1*(V2*(2*(x2*x2) - 2*x1*x5 + 2*x6*(-x3 + x4 + 2*x6) - x2*(x3 + x4 + 4*x6)) + 2*VA*(x1*x2 - 2*(x2*x2) - 2*(x5 + x6)*(x4 + 2*x6) + x2*(x3 + 3*x4 + 2*x5 + 6*x6))) - ml2*ml2*(2*VA*(x1*x2 - 2*(x2*x2) - 2*(x3 + 2*x5)*(x5 + x6) + x2*(3*x3 + x4 + 6*x5 + 2*x6)) + V2*(-2*(x2*x2) + x2*(x3 + x4 + 4*x5) + 2*(-(x3*x5) + x4*x5 - 2*(x5*x5) + x1*x6)))) + 2*(Enu*Enu)*HH1*(Mn*Mn)*x1*(-2*V2*x1*(x2*x2*x2*x2) + 4*V2*(x2*x2*x2*x2*x2) + V2*x1*(x2*x2*x2)*x3 + 2*VA*x1*(x2*x2*x2)*x3 - 6*V2*(x2*x2*x2*x2)*x3 - 4*VA*(x2*x2*x2*x2)*x3 + 2*V2*(x2*x2*x2)*(x3*x3) + 4*VA*(x2*x2*x2)*(x3*x3) + V2*x1*(x2*x2*x2)*x4 - 2*VA*x1*(x2*x2*x2)*x4 - 6*V2*(x2*x2*x2*x2)*x4 + 4*VA*(x2*x2*x2*x2)*x4 + 4*V2*(x2*x2*x2)*x3*x4 + 2*V2*(x2*x2*x2)*(x4*x4) - 4*VA*(x2*x2*x2)*(x4*x4) + 4*V2*(x1*x1)*(x2*x2)*x5 + 8*VA*(x1*x1)*(x2*x2)*x5 - 6*V2*x1*(x2*x2*x2)*x5 - 12*VA*x1*(x2*x2*x2)*x5 - 8*V2*(x2*x2*x2*x2)*x5 - 8*VA*(x2*x2*x2*x2)*x5 + 8*V2*x1*(x2*x2)*x3*x5 + 16*VA*x1*(x2*x2)*x3*x5 + 6*V2*(x2*x2*x2)*x3*x5 + 12*VA*(x2*x2*x2)*x3*x5 + 8*V2*x1*(x2*x2)*x4*x5 + 16*VA*x1*(x2*x2)*x4*x5 + 6*V2*(x2*x2*x2)*x4*x5 - 4*VA*(x2*x2*x2)*x4*x5 + 2*V2*(x2*x2)*x3*x4*x5 + 4*VA*(x2*x2)*x3*x4*x5 - 2*V2*(x2*x2)*(x4*x4)*x5 + 4*VA*(x2*x2)*(x4*x4)*x5 - 3*V2*(x1*x1)*x2*(x5*x5) - 6*VA*(x1*x1)*x2*(x5*x5) + 12*V2*x1*(x2*x2)*(x5*x5) + 24*VA*x1*(x2*x2)*(x5*x5) + 8*V2*(x2*x2*x2)*(x5*x5) + 16*VA*(x2*x2*x2)*(x5*x5) - 6*V2*x1*x2*x3*(x5*x5) - 12*VA*x1*x2*x3*(x5*x5) - 4*V2*(x2*x2)*x3*(x5*x5) - 8*VA*(x2*x2)*x3*(x5*x5) - 12*V2*x1*x2*x4*(x5*x5) - 24*VA*x1*x2*x4*(x5*x5) - 6*V2*x1*x2*(x5*x5*x5) - 12*VA*x1*x2*(x5*x5*x5) - 4*V2*(x2*x2)*(x5*x5*x5) - 8*VA*(x2*x2)*(x5*x5*x5) + 6*V2*x1*x4*(x5*x5*x5) + 12*VA*x1*x4*(x5*x5*x5) + 4*V2*(x1*x1)*(x2*x2)*x6 - 8*VA*(x1*x1)*(x2*x2)*x6 - 6*V2*x1*(x2*x2*x2)*x6 + 12*VA*x1*(x2*x2*x2)*x6 - 8*V2*(x2*x2*x2*x2)*x6 + 8*VA*(x2*x2*x2*x2)*x6 + 8*V2*x1*(x2*x2)*x3*x6 - 16*VA*x1*(x2*x2)*x3*x6 + 6*V2*(x2*x2*x2)*x3*x6 + 4*VA*(x2*x2*x2)*x3*x6 - 2*V2*(x2*x2)*(x3*x3)*x6 - 4*VA*(x2*x2)*(x3*x3)*x6 + 8*V2*x1*(x2*x2)*x4*x6 - 16*VA*x1*(x2*x2)*x4*x6 + 6*V2*(x2*x2*x2)*x4*x6 - 12*VA*(x2*x2*x2)*x4*x6 + 2*V2*(x2*x2)*x3*x4*x6 - 4*VA*(x2*x2)*x3*x4*x6 + 16*V2*x1*(x2*x2)*x5*x6 + 8*V2*(x2*x2*x2)*x5*x6 - 6*V2*x1*x2*x3*x5*x6 - 12*VA*x1*x2*x3*x5*x6 - 4*V2*(x2*x2)*x3*x5*x6 - 8*VA*(x2*x2)*x3*x5*x6 - 6*V2*x1*x2*x4*x5*x6 + 12*VA*x1*x2*x4*x5*x6 - 4*V2*(x2*x2)*x4*x5*x6 + 8*VA*(x2*x2)*x4*x5*x6 - 18*V2*x1*x2*(x5*x5)*x6 - 36*VA*x1*x2*(x5*x5)*x6 - 4*V2*(x2*x2)*(x5*x5)*x6 - 8*VA*(x2*x2)*(x5*x5)*x6 + 6*V2*x1*x3*(x5*x5)*x6 + 12*VA*x1*x3*(x5*x5)*x6 + 12*V2*x1*(x5*x5*x5)*x6 + 24*VA*x1*(x5*x5*x5)*x6 - 3*V2*(x1*x1)*x2*(x6*x6) + 6*VA*(x1*x1)*x2*(x6*x6) + 12*V2*x1*(x2*x2)*(x6*x6) - 24*VA*x1*(x2*x2)*(x6*x6) + 8*V2*(x2*x2*x2)*(x6*x6) - 16*VA*(x2*x2*x2)*(x6*x6) - 12*V2*x1*x2*x3*(x6*x6) + 24*VA*x1*x2*x3*(x6*x6) - 6*V2*x1*x2*x4*(x6*x6) + 12*VA*x1*x2*x4*(x6*x6) - 4*V2*(x2*x2)*x4*(x6*x6) + 8*VA*(x2*x2)*x4*(x6*x6) - 18*V2*x1*x2*x5*(x6*x6) + 36*VA*x1*x2*x5*(x6*x6) - 4*V2*(x2*x2)*x5*(x6*x6) + 8*VA*(x2*x2)*x5*(x6*x6) + 6*V2*x1*x4*x5*(x6*x6) - 12*VA*x1*x4*x5*(x6*x6) - 6*V2*x1*x2*(x6*x6*x6) + 12*VA*x1*x2*(x6*x6*x6) - 4*V2*(x2*x2)*(x6*x6*x6) + 8*VA*(x2*x2)*(x6*x6*x6) + 6*V2*x1*x3*(x6*x6*x6) - 12*VA*x1*x3*(x6*x6*x6) + 12*V2*x1*x5*(x6*x6*x6) - 24*VA*x1*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1)*ml2*V2*(x2*x2)*(-x2 + x5 + x6) + ml2*ml2*ml2*ml2*(x2*x2)*(2*VA*(x2 - x5 - x6) + V2*(-x5 + x6)) + ml1*ml1*ml1*ml1*(x2*x2)*(V2*(x5 - x6) + 2*VA*(-x2 + x5 + x6)) - 2*ml1*ml2*V2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6)))) - ml2*ml2*(V2*(-2*(x2*x2*x2*x2) + x2*x2*x2*(x3 + x4 + 4*x5) + x2*x2*(-2*x5*(x3 - x4 + 2*x5) + x1*(3*x5 - x6)) + 6*x1*x5*(x5 - x6)*x6 + 3*x1*x2*(-(x5*x5) + x6*x6)) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x3 + 2*x5)*(x5 + x6) + x2*(3*x3 + x4 + 6*x5 + 2*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)))) + ml1*ml1*(V2*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x5 + x6)*(x4 + 2*x6) + x2*(x3 + 3*x4 + 2*x5 + 6*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6))))) - 2*Enu*HH1*Mn*x1*x2*(-2*V2*x1*(x2*x2*x2*x2) + 4*V2*(x2*x2*x2*x2*x2) + V2*x1*(x2*x2*x2)*x3 + 2*VA*x1*(x2*x2*x2)*x3 - 6*V2*(x2*x2*x2*x2)*x3 - 4*VA*(x2*x2*x2*x2)*x3 + 2*V2*(x2*x2*x2)*(x3*x3) + 4*VA*(x2*x2*x2)*(x3*x3) + V2*x1*(x2*x2*x2)*x4 - 2*VA*x1*(x2*x2*x2)*x4 - 6*V2*(x2*x2*x2*x2)*x4 + 4*VA*(x2*x2*x2*x2)*x4 + 4*V2*(x2*x2*x2)*x3*x4 + 2*V2*(x2*x2*x2)*(x4*x4) - 4*VA*(x2*x2*x2)*(x4*x4) + 4*V2*(x1*x1)*(x2*x2)*x5 + 8*VA*(x1*x1)*(x2*x2)*x5 - 6*V2*x1*(x2*x2*x2)*x5 - 12*VA*x1*(x2*x2*x2)*x5 - 8*V2*(x2*x2*x2*x2)*x5 - 8*VA*(x2*x2*x2*x2)*x5 + 8*V2*x1*(x2*x2)*x3*x5 + 16*VA*x1*(x2*x2)*x3*x5 + 6*V2*(x2*x2*x2)*x3*x5 + 12*VA*(x2*x2*x2)*x3*x5 + 8*V2*x1*(x2*x2)*x4*x5 + 16*VA*x1*(x2*x2)*x4*x5 + 6*V2*(x2*x2*x2)*x4*x5 - 4*VA*(x2*x2*x2)*x4*x5 + 2*V2*(x2*x2)*x3*x4*x5 + 4*VA*(x2*x2)*x3*x4*x5 - 2*V2*(x2*x2)*(x4*x4)*x5 + 4*VA*(x2*x2)*(x4*x4)*x5 - 3*V2*(x1*x1)*x2*(x5*x5) - 6*VA*(x1*x1)*x2*(x5*x5) + 12*V2*x1*(x2*x2)*(x5*x5) + 24*VA*x1*(x2*x2)*(x5*x5) + 8*V2*(x2*x2*x2)*(x5*x5) + 16*VA*(x2*x2*x2)*(x5*x5) - 6*V2*x1*x2*x3*(x5*x5) - 12*VA*x1*x2*x3*(x5*x5) - 4*V2*(x2*x2)*x3*(x5*x5) - 8*VA*(x2*x2)*x3*(x5*x5) - 12*V2*x1*x2*x4*(x5*x5) - 24*VA*x1*x2*x4*(x5*x5) - 6*V2*x1*x2*(x5*x5*x5) - 12*VA*x1*x2*(x5*x5*x5) - 4*V2*(x2*x2)*(x5*x5*x5) - 8*VA*(x2*x2)*(x5*x5*x5) + 6*V2*x1*x4*(x5*x5*x5) + 12*VA*x1*x4*(x5*x5*x5) + 4*V2*(x1*x1)*(x2*x2)*x6 - 8*VA*(x1*x1)*(x2*x2)*x6 - 6*V2*x1*(x2*x2*x2)*x6 + 12*VA*x1*(x2*x2*x2)*x6 - 8*V2*(x2*x2*x2*x2)*x6 + 8*VA*(x2*x2*x2*x2)*x6 + 8*V2*x1*(x2*x2)*x3*x6 - 16*VA*x1*(x2*x2)*x3*x6 + 6*V2*(x2*x2*x2)*x3*x6 + 4*VA*(x2*x2*x2)*x3*x6 - 2*V2*(x2*x2)*(x3*x3)*x6 - 4*VA*(x2*x2)*(x3*x3)*x6 + 8*V2*x1*(x2*x2)*x4*x6 - 16*VA*x1*(x2*x2)*x4*x6 + 6*V2*(x2*x2*x2)*x4*x6 - 12*VA*(x2*x2*x2)*x4*x6 + 2*V2*(x2*x2)*x3*x4*x6 - 4*VA*(x2*x2)*x3*x4*x6 + 16*V2*x1*(x2*x2)*x5*x6 + 8*V2*(x2*x2*x2)*x5*x6 - 6*V2*x1*x2*x3*x5*x6 - 12*VA*x1*x2*x3*x5*x6 - 4*V2*(x2*x2)*x3*x5*x6 - 8*VA*(x2*x2)*x3*x5*x6 - 6*V2*x1*x2*x4*x5*x6 + 12*VA*x1*x2*x4*x5*x6 - 4*V2*(x2*x2)*x4*x5*x6 + 8*VA*(x2*x2)*x4*x5*x6 - 18*V2*x1*x2*(x5*x5)*x6 - 36*VA*x1*x2*(x5*x5)*x6 - 4*V2*(x2*x2)*(x5*x5)*x6 - 8*VA*(x2*x2)*(x5*x5)*x6 + 6*V2*x1*x3*(x5*x5)*x6 + 12*VA*x1*x3*(x5*x5)*x6 + 12*V2*x1*(x5*x5*x5)*x6 + 24*VA*x1*(x5*x5*x5)*x6 - 3*V2*(x1*x1)*x2*(x6*x6) + 6*VA*(x1*x1)*x2*(x6*x6) + 12*V2*x1*(x2*x2)*(x6*x6) - 24*VA*x1*(x2*x2)*(x6*x6) + 8*V2*(x2*x2*x2)*(x6*x6) - 16*VA*(x2*x2*x2)*(x6*x6) - 12*V2*x1*x2*x3*(x6*x6) + 24*VA*x1*x2*x3*(x6*x6) - 6*V2*x1*x2*x4*(x6*x6) + 12*VA*x1*x2*x4*(x6*x6) - 4*V2*(x2*x2)*x4*(x6*x6) + 8*VA*(x2*x2)*x4*(x6*x6) - 18*V2*x1*x2*x5*(x6*x6) + 36*VA*x1*x2*x5*(x6*x6) - 4*V2*(x2*x2)*x5*(x6*x6) + 8*VA*(x2*x2)*x5*(x6*x6) + 6*V2*x1*x4*x5*(x6*x6) - 12*VA*x1*x4*x5*(x6*x6) - 6*V2*x1*x2*(x6*x6*x6) + 12*VA*x1*x2*(x6*x6*x6) - 4*V2*(x2*x2)*(x6*x6*x6) + 8*VA*(x2*x2)*(x6*x6*x6) + 6*V2*x1*x3*(x6*x6*x6) - 12*VA*x1*x3*(x6*x6*x6) + 12*V2*x1*x5*(x6*x6*x6) - 24*VA*x1*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1)*ml2*V2*(x2*x2)*(-x2 + x5 + x6) + ml2*ml2*ml2*ml2*(x2*x2)*(2*VA*(x2 - x5 - x6) + V2*(-x5 + x6)) + ml1*ml1*ml1*ml1*(x2*x2)*(V2*(x5 - x6) + 2*VA*(-x2 + x5 + x6)) - 2*ml1*ml2*V2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6)))) - ml2*ml2*(V2*(-2*(x2*x2*x2*x2) + x2*x2*x2*(x3 + x4 + 4*x5) + x2*x2*(-2*x5*(x3 - x4 + 2*x5) + x1*(3*x5 - x6)) + 6*x1*x5*(x5 - x6)*x6 + 3*x1*x2*(-(x5*x5) + x6*x6)) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x3 + 2*x5)*(x5 + x6) + x2*(3*x3 + x4 + 6*x5 + 2*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)))) + ml1*ml1*(V2*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x5 + x6)*(x4 + 2*x6) + x2*(x3 + 3*x4 + 2*x5 + 6*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6))))) + HH1*(x2*x2)*(4*(Mn*Mn)*V2*x1*(x2*x2*x2*x2) - 8*(Mn*Mn)*V2*(x2*x2*x2*x2*x2) - 2*(Mn*Mn)*V2*x1*(x2*x2*x2)*x3 - 4*(Mn*Mn)*VA*x1*(x2*x2*x2)*x3 + 12*(Mn*Mn)*V2*(x2*x2*x2*x2)*x3 + 8*(Mn*Mn)*VA*(x2*x2*x2*x2)*x3 - 4*(Mn*Mn)*V2*(x2*x2*x2)*(x3*x3) - 8*(Mn*Mn)*VA*(x2*x2*x2)*(x3*x3) - 2*(Mn*Mn)*V2*x1*(x2*x2*x2)*x4 + 4*(Mn*Mn)*VA*x1*(x2*x2*x2)*x4 + 12*(Mn*Mn)*V2*(x2*x2*x2*x2)*x4 - 8*(Mn*Mn)*VA*(x2*x2*x2*x2)*x4 - 8*(Mn*Mn)*V2*(x2*x2*x2)*x3*x4 - 4*(Mn*Mn)*V2*(x2*x2*x2)*(x4*x4) + 8*(Mn*Mn)*VA*(x2*x2*x2)*(x4*x4) - 4*(Mn*Mn)*V2*(x1*x1)*(x2*x2)*x5 - 8*(Mn*Mn)*VA*(x1*x1)*(x2*x2)*x5 + V2*(x1*x1*x1)*(x2*x2)*x5 + 2*VA*(x1*x1*x1)*(x2*x2)*x5 + 4*(Mn*Mn)*V2*x1*(x2*x2*x2)*x5 + 8*(Mn*Mn)*VA*x1*(x2*x2*x2)*x5 - 2*V2*(x1*x1)*(x2*x2*x2)*x5 - 4*VA*(x1*x1)*(x2*x2*x2)*x5 + 16*(Mn*Mn)*V2*(x2*x2*x2*x2)*x5 + 16*(Mn*Mn)*VA*(x2*x2*x2*x2)*x5 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x3*x5 - 16*(Mn*Mn)*VA*x1*(x2*x2)*x3*x5 + 2*V2*(x1*x1)*(x2*x2)*x3*x5 + 4*VA*(x1*x1)*(x2*x2)*x3*x5 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x3*x5 - 24*(Mn*Mn)*VA*(x2*x2*x2)*x3*x5 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x4*x5 - 16*(Mn*Mn)*VA*x1*(x2*x2)*x4*x5 + 2*V2*(x1*x1)*(x2*x2)*x4*x5 + 4*VA*(x1*x1)*(x2*x2)*x4*x5 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x4*x5 + 8*(Mn*Mn)*VA*(x2*x2*x2)*x4*x5 - 4*(Mn*Mn)*V2*(x2*x2)*x3*x4*x5 - 8*(Mn*Mn)*VA*(x2*x2)*x3*x4*x5 + 4*(Mn*Mn)*V2*(x2*x2)*(x4*x4)*x5 - 8*(Mn*Mn)*VA*(x2*x2)*(x4*x4)*x5 + 2*(Mn*Mn)*V2*(x1*x1)*x2*(x5*x5) + 4*(Mn*Mn)*VA*(x1*x1)*x2*(x5*x5) - V2*(x1*x1*x1)*x2*(x5*x5) - 2*VA*(x1*x1*x1)*x2*(x5*x5) - 8*(Mn*Mn)*V2*x1*(x2*x2)*(x5*x5) - 16*(Mn*Mn)*VA*x1*(x2*x2)*(x5*x5) + 4*V2*(x1*x1)*(x2*x2)*(x5*x5) + 8*VA*(x1*x1)*(x2*x2)*(x5*x5) - 16*(Mn*Mn)*V2*(x2*x2*x2)*(x5*x5) - 32*(Mn*Mn)*VA*(x2*x2*x2)*(x5*x5) + 4*(Mn*Mn)*V2*x1*x2*x3*(x5*x5) + 8*(Mn*Mn)*VA*x1*x2*x3*(x5*x5) - 2*V2*(x1*x1)*x2*x3*(x5*x5) - 4*VA*(x1*x1)*x2*x3*(x5*x5) + 8*(Mn*Mn)*V2*(x2*x2)*x3*(x5*x5) + 16*(Mn*Mn)*VA*(x2*x2)*x3*(x5*x5) + 8*(Mn*Mn)*V2*x1*x2*x4*(x5*x5) + 16*(Mn*Mn)*VA*x1*x2*x4*(x5*x5) - 4*V2*(x1*x1)*x2*x4*(x5*x5) - 8*VA*(x1*x1)*x2*x4*(x5*x5) + 4*(Mn*Mn)*V2*x1*x2*(x5*x5*x5) + 8*(Mn*Mn)*VA*x1*x2*(x5*x5*x5) - 2*V2*(x1*x1)*x2*(x5*x5*x5) - 4*VA*(x1*x1)*x2*(x5*x5*x5) + 8*(Mn*Mn)*V2*(x2*x2)*(x5*x5*x5) + 16*(Mn*Mn)*VA*(x2*x2)*(x5*x5*x5) - 4*(Mn*Mn)*V2*x1*x4*(x5*x5*x5) - 8*(Mn*Mn)*VA*x1*x4*(x5*x5*x5) + 2*V2*(x1*x1)*x4*(x5*x5*x5) + 4*VA*(x1*x1)*x4*(x5*x5*x5) + 4*(ml1*ml1*ml1)*ml2*(Mn*Mn)*V2*(x2*x2)*(x2 - x5 - x6) - 4*(Mn*Mn)*V2*(x1*x1)*(x2*x2)*x6 + 8*(Mn*Mn)*VA*(x1*x1)*(x2*x2)*x6 + V2*(x1*x1*x1)*(x2*x2)*x6 - 2*VA*(x1*x1*x1)*(x2*x2)*x6 + 4*(Mn*Mn)*V2*x1*(x2*x2*x2)*x6 - 8*(Mn*Mn)*VA*x1*(x2*x2*x2)*x6 - 2*V2*(x1*x1)*(x2*x2*x2)*x6 + 4*VA*(x1*x1)*(x2*x2*x2)*x6 + 16*(Mn*Mn)*V2*(x2*x2*x2*x2)*x6 - 16*(Mn*Mn)*VA*(x2*x2*x2*x2)*x6 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x3*x6 + 16*(Mn*Mn)*VA*x1*(x2*x2)*x3*x6 + 2*V2*(x1*x1)*(x2*x2)*x3*x6 - 4*VA*(x1*x1)*(x2*x2)*x3*x6 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x3*x6 - 8*(Mn*Mn)*VA*(x2*x2*x2)*x3*x6 + 4*(Mn*Mn)*V2*(x2*x2)*(x3*x3)*x6 + 8*(Mn*Mn)*VA*(x2*x2)*(x3*x3)*x6 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x4*x6 + 16*(Mn*Mn)*VA*x1*(x2*x2)*x4*x6 + 2*V2*(x1*x1)*(x2*x2)*x4*x6 - 4*VA*(x1*x1)*(x2*x2)*x4*x6 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x4*x6 + 24*(Mn*Mn)*VA*(x2*x2*x2)*x4*x6 - 4*(Mn*Mn)*V2*(x2*x2)*x3*x4*x6 + 8*(Mn*Mn)*VA*(x2*x2)*x3*x4*x6 - 16*(Mn*Mn)*V2*x1*(x2*x2)*x5*x6 + 4*V2*(x1*x1)*(x2*x2)*x5*x6 - 16*(Mn*Mn)*V2*(x2*x2*x2)*x5*x6 + 4*(Mn*Mn)*V2*x1*x2*x3*x5*x6 + 8*(Mn*Mn)*VA*x1*x2*x3*x5*x6 - 2*V2*(x1*x1)*x2*x3*x5*x6 - 4*VA*(x1*x1)*x2*x3*x5*x6 + 8*(Mn*Mn)*V2*(x2*x2)*x3*x5*x6 + 16*(Mn*Mn)*VA*(x2*x2)*x3*x5*x6 + 4*(Mn*Mn)*V2*x1*x2*x4*x5*x6 - 8*(Mn*Mn)*VA*x1*x2*x4*x5*x6 - 2*V2*(x1*x1)*x2*x4*x5*x6 + 4*VA*(x1*x1)*x2*x4*x5*x6 + 8*(Mn*Mn)*V2*(x2*x2)*x4*x5*x6 - 16*(Mn*Mn)*VA*(x2*x2)*x4*x5*x6 + 12*(Mn*Mn)*V2*x1*x2*(x5*x5)*x6 + 24*(Mn*Mn)*VA*x1*x2*(x5*x5)*x6 - 6*V2*(x1*x1)*x2*(x5*x5)*x6 - 12*VA*(x1*x1)*x2*(x5*x5)*x6 + 8*(Mn*Mn)*V2*(x2*x2)*(x5*x5)*x6 + 16*(Mn*Mn)*VA*(x2*x2)*(x5*x5)*x6 - 4*(Mn*Mn)*V2*x1*x3*(x5*x5)*x6 - 8*(Mn*Mn)*VA*x1*x3*(x5*x5)*x6 + 2*V2*(x1*x1)*x3*(x5*x5)*x6 + 4*VA*(x1*x1)*x3*(x5*x5)*x6 - 8*(Mn*Mn)*V2*x1*(x5*x5*x5)*x6 - 16*(Mn*Mn)*VA*x1*(x5*x5*x5)*x6 + 4*V2*(x1*x1)*(x5*x5*x5)*x6 + 8*VA*(x1*x1)*(x5*x5*x5)*x6 + 2*(Mn*Mn)*V2*(x1*x1)*x2*(x6*x6) - 4*(Mn*Mn)*VA*(x1*x1)*x2*(x6*x6) - V2*(x1*x1*x1)*x2*(x6*x6) + 2*VA*(x1*x1*x1)*x2*(x6*x6) - 8*(Mn*Mn)*V2*x1*(x2*x2)*(x6*x6) + 16*(Mn*Mn)*VA*x1*(x2*x2)*(x6*x6) + 4*V2*(x1*x1)*(x2*x2)*(x6*x6) - 8*VA*(x1*x1)*(x2*x2)*(x6*x6) - 16*(Mn*Mn)*V2*(x2*x2*x2)*(x6*x6) + 32*(Mn*Mn)*VA*(x2*x2*x2)*(x6*x6) + 8*(Mn*Mn)*V2*x1*x2*x3*(x6*x6) - 16*(Mn*Mn)*VA*x1*x2*x3*(x6*x6) - 4*V2*(x1*x1)*x2*x3*(x6*x6) + 8*VA*(x1*x1)*x2*x3*(x6*x6) + 4*(Mn*Mn)*V2*x1*x2*x4*(x6*x6) - 8*(Mn*Mn)*VA*x1*x2*x4*(x6*x6) - 2*V2*(x1*x1)*x2*x4*(x6*x6) + 4*VA*(x1*x1)*x2*x4*(x6*x6) + 8*(Mn*Mn)*V2*(x2*x2)*x4*(x6*x6) - 16*(Mn*Mn)*VA*(x2*x2)*x4*(x6*x6) + 12*(Mn*Mn)*V2*x1*x2*x5*(x6*x6) - 24*(Mn*Mn)*VA*x1*x2*x5*(x6*x6) - 6*V2*(x1*x1)*x2*x5*(x6*x6) + 12*VA*(x1*x1)*x2*x5*(x6*x6) + 8*(Mn*Mn)*V2*(x2*x2)*x5*(x6*x6) - 16*(Mn*Mn)*VA*(x2*x2)*x5*(x6*x6) - 4*(Mn*Mn)*V2*x1*x4*x5*(x6*x6) + 8*(Mn*Mn)*VA*x1*x4*x5*(x6*x6) + 2*V2*(x1*x1)*x4*x5*(x6*x6) - 4*VA*(x1*x1)*x4*x5*(x6*x6) + 4*(Mn*Mn)*V2*x1*x2*(x6*x6*x6) - 8*(Mn*Mn)*VA*x1*x2*(x6*x6*x6) - 2*V2*(x1*x1)*x2*(x6*x6*x6) + 4*VA*(x1*x1)*x2*(x6*x6*x6) + 8*(Mn*Mn)*V2*(x2*x2)*(x6*x6*x6) - 16*(Mn*Mn)*VA*(x2*x2)*(x6*x6*x6) - 4*(Mn*Mn)*V2*x1*x3*(x6*x6*x6) + 8*(Mn*Mn)*VA*x1*x3*(x6*x6*x6) + 2*V2*(x1*x1)*x3*(x6*x6*x6) - 4*VA*(x1*x1)*x3*(x6*x6*x6) - 8*(Mn*Mn)*V2*x1*x5*(x6*x6*x6) + 16*(Mn*Mn)*VA*x1*x5*(x6*x6*x6) + 4*V2*(x1*x1)*x5*(x6*x6*x6) - 8*VA*(x1*x1)*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1*ml1)*(Mn*Mn)*(x2*x2)*(2*VA*(x2 - x5 - x6) + V2*(-x5 + x6)) + 2*(ml2*ml2*ml2*ml2)*(Mn*Mn)*(x2*x2)*(V2*(x5 - x6) + 2*VA*(-x2 + x5 + x6)) + 2*ml1*ml2*V2*(2*(ml2*ml2)*(Mn*Mn)*(x2*x2)*(x2 - x5 - x6) - x1*x1*(x2*x2*x2 - 2*(x2*x2)*(x5 + x6) - 2*x5*x6*(x5 + x6) + x2*(x5*x5 + 4*x5*x6 + x6*x6)) + 2*(Mn*Mn)*(-(x1*(x2*x2*x2 + 2*(x2*x2)*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6))))) + ml2*ml2*(x1*x1*(-(V2*(x5 - x6)*(x2*x2 + 2*x5*x6 - x2*(x5 + x6))) - 2*VA*(x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) + 2*(Mn*Mn)*(V2*(-2*(x2*x2*x2*x2) + x2*x2*x2*(x3 + x4 + 4*x5) + 2*x1*x5*(x5 - x6)*x6 + x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(-2*x5*(x3 - x4 + 2*x5) + x1*(x5 + x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x3 + 2*x5)*(x5 + x6) + x2*(3*x3 + x4 + 6*x5 + 2*x6)) + x1*(x2*x2*x2 + x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))))) + ml1*ml1*(x1*x1*(V2*(x5 - x6)*(x2*x2 + 2*x5*x6 - x2*(x5 + x6)) + 2*VA*(x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) - 2*(Mn*Mn)*(V2*(2*(x2*x2*x2*x2) + 2*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + x1*x2*(-(x5*x5) + x6*x6) - x2*x2*(2*(x3 - x4 - 2*x6)*x6 + x1*(x5 + x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x5 + x6)*(x4 + 2*x6) + x2*(x3 + 3*x4 + 2*x5 + 6*x6)) + x1*(x2*x2*x2 + x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6)))))) + A2*(HH2*x1*(x2*x2*x2*x2)*(-2*x1*(x2*x2) + 4*(x2*x2*x2) + x1*x2*x3 - 6*(x2*x2)*x3 + 2*x2*(x3*x3) + x1*x2*x4 - 6*(x2*x2)*x4 + 4*x2*x3*x4 + 2*x2*(x4*x4) + x1*x1*x5 - 8*(x2*x2)*x5 + 2*x1*x3*x5 + 6*x2*x3*x5 + 2*x1*x4*x5 + 6*x2*x4*x5 + 2*x3*x4*x5 - 2*(x4*x4)*x5 + 8*x2*(x5*x5) - 4*x3*(x5*x5) - 4*(x5*x5*x5) + 2*(ml1*ml1*ml1)*ml2*(x2 - x5 - x6) + ml1*ml1*ml1*ml1*(x5 - x6) + x1*x1*x6 - 8*(x2*x2)*x6 + 2*x1*x3*x6 + 6*x2*x3*x6 - 2*(x3*x3)*x6 + 2*x1*x4*x6 + 6*x2*x4*x6 + 2*x3*x4*x6 + 4*x1*x5*x6 + 8*x2*x5*x6 - 4*x3*x5*x6 - 4*x4*x5*x6 - 4*(x5*x5)*x6 + 8*x2*(x6*x6) - 4*x4*(x6*x6) - 4*x5*(x6*x6) - 4*(x6*x6*x6) + ml2*ml2*ml2*ml2*(-x5 + x6) + ml2*ml2*(2*(x2*x2) + 2*x3*x5 - 2*x4*x5 + 4*(x5*x5) - x2*(x3 + x4 + 4*x5) - 2*x1*x6) + 2*ml1*ml2*(-2*x1*x2 + 4*(x2*x2) - 3*x2*x3 - 3*x2*x4 - 4*x2*x5 + x3*x5 + x4*x5 + 2*(x5*x5) + ml2*ml2*(x2 - x5 - x6) - 4*x2*x6 + x3*x6 + x4*x6 + 4*x5*x6 + 2*(x6*x6)) + ml1*ml1*(2*(x2*x2) - 2*x1*x5 + 2*x6*(-x3 + x4 + 2*x6) - x2*(x3 + x4 + 4*x6))) + 2*(Enu*Enu)*HH1*(Mn*Mn)*x1*(-2*x1*(x2*x2*x2*x2) + 4*(x2*x2*x2*x2*x2) + x1*(x2*x2*x2)*x3 - 6*(x2*x2*x2*x2)*x3 + 2*(x2*x2*x2)*(x3*x3) + x1*(x2*x2*x2)*x4 - 6*(x2*x2*x2*x2)*x4 + 4*(x2*x2*x2)*x3*x4 + 2*(x2*x2*x2)*(x4*x4) + 4*(x1*x1)*(x2*x2)*x5 - 6*x1*(x2*x2*x2)*x5 - 8*(x2*x2*x2*x2)*x5 + 8*x1*(x2*x2)*x3*x5 + 6*(x2*x2*x2)*x3*x5 + 8*x1*(x2*x2)*x4*x5 + 6*(x2*x2*x2)*x4*x5 + 2*(x2*x2)*x3*x4*x5 - 2*(x2*x2)*(x4*x4)*x5 - 3*(x1*x1)*x2*(x5*x5) + 12*x1*(x2*x2)*(x5*x5) + 8*(x2*x2*x2)*(x5*x5) - 6*x1*x2*x3*(x5*x5) - 4*(x2*x2)*x3*(x5*x5) - 12*x1*x2*x4*(x5*x5) - 6*x1*x2*(x5*x5*x5) - 4*(x2*x2)*(x5*x5*x5) + 6*x1*x4*(x5*x5*x5) + 2*(ml1*ml1*ml1)*ml2*(x2*x2)*(x2 - x5 - x6) + ml1*ml1*ml1*ml1*(x2*x2)*(x5 - x6) + 4*(x1*x1)*(x2*x2)*x6 - 6*x1*(x2*x2*x2)*x6 - 8*(x2*x2*x2*x2)*x6 + 8*x1*(x2*x2)*x3*x6 + 6*(x2*x2*x2)*x3*x6 - 2*(x2*x2)*(x3*x3)*x6 + 8*x1*(x2*x2)*x4*x6 + 6*(x2*x2*x2)*x4*x6 + 2*(x2*x2)*x3*x4*x6 + 16*x1*(x2*x2)*x5*x6 + 8*(x2*x2*x2)*x5*x6 - 6*x1*x2*x3*x5*x6 - 4*(x2*x2)*x3*x5*x6 - 6*x1*x2*x4*x5*x6 - 4*(x2*x2)*x4*x5*x6 - 18*x1*x2*(x5*x5)*x6 - 4*(x2*x2)*(x5*x5)*x6 + 6*x1*x3*(x5*x5)*x6 + 12*x1*(x5*x5*x5)*x6 - 3*(x1*x1)*x2*(x6*x6) + 12*x1*(x2*x2)*(x6*x6) + 8*(x2*x2*x2)*(x6*x6) - 12*x1*x2*x3*(x6*x6) - 6*x1*x2*x4*(x6*x6) - 4*(x2*x2)*x4*(x6*x6) - 18*x1*x2*x5*(x6*x6) - 4*(x2*x2)*x5*(x6*x6) + 6*x1*x4*x5*(x6*x6) - 6*x1*x2*(x6*x6*x6) - 4*(x2*x2)*(x6*x6*x6) + 6*x1*x3*(x6*x6*x6) + 12*x1*x5*(x6*x6*x6) + ml2*ml2*ml2*ml2*(x2*x2)*(-x5 + x6) + ml2*ml2*(2*(x2*x2*x2*x2) - x2*x2*x2*(x3 + x4 + 4*x5) + 6*x1*x5*x6*(-x5 + x6) + 3*x1*x2*(x5*x5 - x6*x6) + x2*x2*(2*x5*(x3 - x4 + 2*x5) + x1*(-3*x5 + x6))) + ml1*ml1*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*ml1*ml2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6))))) - 2*Enu*HH1*Mn*x1*x2*(-2*x1*(x2*x2*x2*x2) + 4*(x2*x2*x2*x2*x2) + x1*(x2*x2*x2)*x3 - 6*(x2*x2*x2*x2)*x3 + 2*(x2*x2*x2)*(x3*x3) + x1*(x2*x2*x2)*x4 - 6*(x2*x2*x2*x2)*x4 + 4*(x2*x2*x2)*x3*x4 + 2*(x2*x2*x2)*(x4*x4) + 4*(x1*x1)*(x2*x2)*x5 - 6*x1*(x2*x2*x2)*x5 - 8*(x2*x2*x2*x2)*x5 + 8*x1*(x2*x2)*x3*x5 + 6*(x2*x2*x2)*x3*x5 + 8*x1*(x2*x2)*x4*x5 + 6*(x2*x2*x2)*x4*x5 + 2*(x2*x2)*x3*x4*x5 - 2*(x2*x2)*(x4*x4)*x5 - 3*(x1*x1)*x2*(x5*x5) + 12*x1*(x2*x2)*(x5*x5) + 8*(x2*x2*x2)*(x5*x5) - 6*x1*x2*x3*(x5*x5) - 4*(x2*x2)*x3*(x5*x5) - 12*x1*x2*x4*(x5*x5) - 6*x1*x2*(x5*x5*x5) - 4*(x2*x2)*(x5*x5*x5) + 6*x1*x4*(x5*x5*x5) + 2*(ml1*ml1*ml1)*ml2*(x2*x2)*(x2 - x5 - x6) + ml1*ml1*ml1*ml1*(x2*x2)*(x5 - x6) + 4*(x1*x1)*(x2*x2)*x6 - 6*x1*(x2*x2*x2)*x6 - 8*(x2*x2*x2*x2)*x6 + 8*x1*(x2*x2)*x3*x6 + 6*(x2*x2*x2)*x3*x6 - 2*(x2*x2)*(x3*x3)*x6 + 8*x1*(x2*x2)*x4*x6 + 6*(x2*x2*x2)*x4*x6 + 2*(x2*x2)*x3*x4*x6 + 16*x1*(x2*x2)*x5*x6 + 8*(x2*x2*x2)*x5*x6 - 6*x1*x2*x3*x5*x6 - 4*(x2*x2)*x3*x5*x6 - 6*x1*x2*x4*x5*x6 - 4*(x2*x2)*x4*x5*x6 - 18*x1*x2*(x5*x5)*x6 - 4*(x2*x2)*(x5*x5)*x6 + 6*x1*x3*(x5*x5)*x6 + 12*x1*(x5*x5*x5)*x6 - 3*(x1*x1)*x2*(x6*x6) + 12*x1*(x2*x2)*(x6*x6) + 8*(x2*x2*x2)*(x6*x6) - 12*x1*x2*x3*(x6*x6) - 6*x1*x2*x4*(x6*x6) - 4*(x2*x2)*x4*(x6*x6) - 18*x1*x2*x5*(x6*x6) - 4*(x2*x2)*x5*(x6*x6) + 6*x1*x4*x5*(x6*x6) - 6*x1*x2*(x6*x6*x6) - 4*(x2*x2)*(x6*x6*x6) + 6*x1*x3*(x6*x6*x6) + 12*x1*x5*(x6*x6*x6) + ml2*ml2*ml2*ml2*(x2*x2)*(-x5 + x6) + ml2*ml2*(2*(x2*x2*x2*x2) - x2*x2*x2*(x3 + x4 + 4*x5) + 6*x1*x5*x6*(-x5 + x6) + 3*x1*x2*(x5*x5 - x6*x6) + x2*x2*(2*x5*(x3 - x4 + 2*x5) + x1*(-3*x5 + x6))) + ml1*ml1*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*ml1*ml2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6))))) + HH1*(x2*x2)*(-4*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2*x2) + 4*(Mn*Mn)*x1*(x2*x2*x2*x2) - 8*(Mn*Mn)*(x2*x2*x2*x2*x2) + 2*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2)*x3 - 2*(Mn*Mn)*x1*(x2*x2*x2)*x3 + 12*(Mn*Mn)*(x2*x2*x2*x2)*x3 - 4*(Mn*Mn)*(x2*x2*x2)*(x3*x3) + 2*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2)*x4 - 2*(Mn*Mn)*x1*(x2*x2*x2)*x4 + 12*(Mn*Mn)*(x2*x2*x2*x2)*x4 - 8*(Mn*Mn)*(x2*x2*x2)*x3*x4 - 4*(Mn*Mn)*(x2*x2*x2)*(x4*x4) + 2*(ml2*ml2*ml2*ml2)*(Mn*Mn)*(x2*x2)*x5 + 2*(ml2*ml2)*(Mn*Mn)*x1*(x2*x2)*x5 - ml2*ml2*(x1*x1)*(x2*x2)*x5 - 4*(Mn*Mn)*(x1*x1)*(x2*x2)*x5 + x1*x1*x1*(x2*x2)*x5 + 8*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2)*x5 + 4*(Mn*Mn)*x1*(x2*x2*x2)*x5 - 2*(x1*x1)*(x2*x2*x2)*x5 + 16*(Mn*Mn)*(x2*x2*x2*x2)*x5 - 4*(ml2*ml2)*(Mn*Mn)*(x2*x2)*x3*x5 - 8*(Mn*Mn)*x1*(x2*x2)*x3*x5 + 2*(x1*x1)*(x2*x2)*x3*x5 - 12*(Mn*Mn)*(x2*x2*x2)*x3*x5 + 4*(ml2*ml2)*(Mn*Mn)*(x2*x2)*x4*x5 - 8*(Mn*Mn)*x1*(x2*x2)*x4*x5 + 2*(x1*x1)*(x2*x2)*x4*x5 - 12*(Mn*Mn)*(x2*x2*x2)*x4*x5 - 4*(Mn*Mn)*(x2*x2)*x3*x4*x5 + 4*(Mn*Mn)*(x2*x2)*(x4*x4)*x5 - 2*(ml2*ml2)*(Mn*Mn)*x1*x2*(x5*x5) + ml2*ml2*(x1*x1)*x2*(x5*x5) + 2*(Mn*Mn)*(x1*x1)*x2*(x5*x5) - x1*x1*x1*x2*(x5*x5) - 8*(ml2*ml2)*(Mn*Mn)*(x2*x2)*(x5*x5) - 8*(Mn*Mn)*x1*(x2*x2)*(x5*x5) + 4*(x1*x1)*(x2*x2)*(x5*x5) - 16*(Mn*Mn)*(x2*x2*x2)*(x5*x5) + 4*(Mn*Mn)*x1*x2*x3*(x5*x5) - 2*(x1*x1)*x2*x3*(x5*x5) + 8*(Mn*Mn)*(x2*x2)*x3*(x5*x5) + 8*(Mn*Mn)*x1*x2*x4*(x5*x5) - 4*(x1*x1)*x2*x4*(x5*x5) + 4*(Mn*Mn)*x1*x2*(x5*x5*x5) - 2*(x1*x1)*x2*(x5*x5*x5) + 8*(Mn*Mn)*(x2*x2)*(x5*x5*x5) - 4*(Mn*Mn)*x1*x4*(x5*x5*x5) + 2*(x1*x1)*x4*(x5*x5*x5) - 2*(ml2*ml2*ml2*ml2)*(Mn*Mn)*(x2*x2)*x6 + 2*(ml2*ml2)*(Mn*Mn)*x1*(x2*x2)*x6 + ml2*ml2*(x1*x1)*(x2*x2)*x6 - 4*(Mn*Mn)*(x1*x1)*(x2*x2)*x6 + x1*x1*x1*(x2*x2)*x6 + 4*(Mn*Mn)*x1*(x2*x2*x2)*x6 - 2*(x1*x1)*(x2*x2*x2)*x6 + 16*(Mn*Mn)*(x2*x2*x2*x2)*x6 - 8*(Mn*Mn)*x1*(x2*x2)*x3*x6 + 2*(x1*x1)*(x2*x2)*x3*x6 - 12*(Mn*Mn)*(x2*x2*x2)*x3*x6 + 4*(Mn*Mn)*(x2*x2)*(x3*x3)*x6 - 8*(Mn*Mn)*x1*(x2*x2)*x4*x6 + 2*(x1*x1)*(x2*x2)*x4*x6 - 12*(Mn*Mn)*(x2*x2*x2)*x4*x6 - 4*(Mn*Mn)*(x2*x2)*x3*x4*x6 - 16*(Mn*Mn)*x1*(x2*x2)*x5*x6 + 4*(x1*x1)*(x2*x2)*x5*x6 - 16*(Mn*Mn)*(x2*x2*x2)*x5*x6 + 4*(Mn*Mn)*x1*x2*x3*x5*x6 - 2*(x1*x1)*x2*x3*x5*x6 + 8*(Mn*Mn)*(x2*x2)*x3*x5*x6 + 4*(Mn*Mn)*x1*x2*x4*x5*x6 - 2*(x1*x1)*x2*x4*x5*x6 + 8*(Mn*Mn)*(x2*x2)*x4*x5*x6 + 4*(ml2*ml2)*(Mn*Mn)*x1*(x5*x5)*x6 - 2*(ml2*ml2)*(x1*x1)*(x5*x5)*x6 + 12*(Mn*Mn)*x1*x2*(x5*x5)*x6 - 6*(x1*x1)*x2*(x5*x5)*x6 + 8*(Mn*Mn)*(x2*x2)*(x5*x5)*x6 - 4*(Mn*Mn)*x1*x3*(x5*x5)*x6 + 2*(x1*x1)*x3*(x5*x5)*x6 - 8*(Mn*Mn)*x1*(x5*x5*x5)*x6 + 4*(x1*x1)*(x5*x5*x5)*x6 + 2*(ml2*ml2)*(Mn*Mn)*x1*x2*(x6*x6) - ml2*ml2*(x1*x1)*x2*(x6*x6) + 2*(Mn*Mn)*(x1*x1)*x2*(x6*x6) - x1*x1*x1*x2*(x6*x6) - 8*(Mn*Mn)*x1*(x2*x2)*(x6*x6) + 4*(x1*x1)*(x2*x2)*(x6*x6) - 16*(Mn*Mn)*(x2*x2*x2)*(x6*x6) + 8*(Mn*Mn)*x1*x2*x3*(x6*x6) - 4*(x1*x1)*x2*x3*(x6*x6) + 4*(Mn*Mn)*x1*x2*x4*(x6*x6) - 2*(x1*x1)*x2*x4*(x6*x6) + 8*(Mn*Mn)*(x2*x2)*x4*(x6*x6) - 4*(ml2*ml2)*(Mn*Mn)*x1*x5*(x6*x6) + 2*(ml2*ml2)*(x1*x1)*x5*(x6*x6) + 12*(Mn*Mn)*x1*x2*x5*(x6*x6) - 6*(x1*x1)*x2*x5*(x6*x6) + 8*(Mn*Mn)*(x2*x2)*x5*(x6*x6) - 4*(Mn*Mn)*x1*x4*x5*(x6*x6) + 2*(x1*x1)*x4*x5*(x6*x6) + 4*(Mn*Mn)*x1*x2*(x6*x6*x6) - 2*(x1*x1)*x2*(x6*x6*x6) + 8*(Mn*Mn)*(x2*x2)*(x6*x6*x6) - 4*(Mn*Mn)*x1*x3*(x6*x6*x6) + 2*(x1*x1)*x3*(x6*x6*x6) - 8*(Mn*Mn)*x1*x5*(x6*x6*x6) + 4*(x1*x1)*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1*ml1)*(Mn*Mn)*(x2*x2)*(-x5 + x6) + 4*(ml1*ml1*ml1)*ml2*(Mn*Mn)*(x2*x2)*(-x2 + x5 + x6) + ml1*ml1*(x1*x1*(x5 - x6)*(x2*x2 + 2*x5*x6 - x2*(x5 + x6)) + 2*(Mn*Mn)*(-2*(x2*x2*x2*x2) + 2*x1*x5*x6*(-x5 + x6) + x2*x2*x2*(x3 + x4 + 4*x6) + x1*x2*(x5*x5 - x6*x6) + x2*x2*(2*(x3 - x4 - 2*x6)*x6 + x1*(x5 + x6)))) - 2*ml1*ml2*(2*(ml2*ml2)*(Mn*Mn)*(x2*x2)*(x2 - x5 - x6) - x1*x1*(x2*x2*x2 - 2*(x2*x2)*(x5 + x6) - 2*x5*x6*(x5 + x6) + x2*(x5*x5 + 4*x5*x6 + x6*x6)) + 2*(Mn*Mn)*(-(x1*(x2*x2*x2 + 2*(x2*x2)*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6))))))))))/(Enu*Enu*(Mn*Mn)*(x1*x1)*(x2*x2*x2*x2)*((x1 + 2*x3)*(x1 + 2*x3))*((x1 + 2*x4)*(x1 + 2*x4)));

  // FERMI BLOCKING
  if(myMC->PAULI_BLOCKING == W_BLOCKING){
    dsigma *= F_PAULI_BLOCK(x1);
  }


  ///////////////////////////////////////////////////////////////////////////////////
  // JACOBIANS AND PREFACTORS 
  dsigma *= Jacob;
  
  // Changing units to zeptobarn = zb = 1e-45 cm^2
  dsigma *= (GeV2_to_cm2*1e45);
  

  myMC->xphys[0] = x1;
  myMC->xphys[1] = x2;
  myMC->xphys[2] = x3;
  myMC->xphys[3] = x4;
  myMC->xphys[4] = x5;
  myMC->xphys[5] = x6;
  myMC->xphys[6] = x7;
  myMC->xphys[7] = x8;
  myMC->xphys[8] = Enu;
  myMC->xphys[9] = dsigma;

  return 0;
}


int DifFromIntToPhysical(const cubareal xx[], void *MC){

  long double u1_s    = xx[0];
  long double u2_s    = xx[1];
  long double u3_s    = xx[2];
  long double PHI2_s  = xx[3];
  long double x5_s    = xx[4];
  long double u6_s    = xx[5];
  long double x7_s    = xx[6];
  long double x8_s    = xx[7];
  long double Enu_s   = xx[8];

  // Random angles needed to specify 4 vectors
  long double x7 = x7_s*2*M_PI;
  long double x8 = x8_s*2*M_PI;

  //////////////////////////////////////////////////
  // get data from user

  tridentMC* myMC = (tridentMC *)MC;

  long double ml1 = myMC->ml1;
  long double ml2 = myMC->ml2;
  long double A = myMC->A;
  long double Z = myMC->Z;
  long double Mn = myMC->Mn;
  long double mzprime = myMC->mzprime;


  long double Diag11 = myMC->terms[0];
  long double Diag22 = myMC->terms[1];
  long double Diag12 = myMC->terms[2];
  long double V2     = myMC->terms[3];
  long double A2     = myMC->terms[4];
  long double VA     = myMC->terms[5];
  long double BSM    = myMC->terms[6];

  long double Enu = Enu_s*(myMC->Emax - myMC->Emin) + myMC->Emin;


  //////////////////////////////////
  long double x1_u = 2*SQR(Enu) / ( 1 + 2 * Enu/Mn ) * ( 1 - SQR(ml1+ml2)/2.0/SQR(Enu) * (1 + Enu/Mn) +
                  sqrt( SQR(1 - SQR(ml1 + ml2)/2.0/SQR(Enu)*(1+Enu/Mn) ) 
                    - pow(ml1 + ml2, 4)/4.0/pow(Enu,4)*(1 + 2*Enu/Mn) )  );
  long double x1_l = pow(ml1+ml2,4)/(x1_u)/(1.0+2*Enu/Mn);

  long double u1_l = log(x1_l); 
  long double u1_u = log(x1_u); 
  
  long double u1 = u1_s*(u1_u-u1_l) + u1_l;

  long double x1 = exp(u1);

  //////////////////////////////////
  long double x2_l = 1.0/2.0 * ( x1 + SQR(ml1 + ml2) );
  long double x2_u = Enu * (sqrt( x1 + SQR(x1)/4.0/SQR(Mn) ) - x1/2.0/Mn );
  // long double x2 = x2_s*(x2_u-x2_l) + x2_l;

  long double u2_u = (1.0+2.0*Enu/Mn)*(x1_u-x1)*(x1-x1_l)  / ( SQR(ml1 + ml2) + x1*(1+Enu/Mn) + 2*Enu*sqrt(x1+SQR(x1)/4.0/SQR(Mn)) );
  long double u2_l = 0.0;

  long double u2 = u2_s*(u2_u-u2_l) + u2_l;

  long double x2 = 0.5*(u2 + SQR(ml1 + ml2) + x1);

  /////////////////
  long double x5_l = x2/2.0/(2*x2-x1) * (  2*x2 - x1 + SQR(ml2) - SQR(ml1)  
        - sqrt( SQR( 2*x2-x1 + SQR(ml2) - SQR(ml1) ) - 4*SQR(ml2)*(2*x2-x1) ) );
  
  long double x5_u = x2/2.0/(2*x2-x1) * (  2*x2 - x1 + SQR(ml2) - SQR(ml1)  
        + sqrt( SQR( 2*x2-x1 + SQR(ml2) - SQR(ml1) ) - 4*SQR(ml2)*(2*x2-x1) ) );
  

  long double x5 = x5_s*(x5_u-x5_l) + x5_l; 

  ////////////////
  long double x3_l = SQR(ml2)/2.0 * x2/x5 - x1/2.0 * x5/x2;

  long double x3_u = 0.5*(2*x2 - x1 + SQR(ml2) - SQR(ml1) )  - x5;
  
  long double u3_l = log(2*x3_l + x1);
  long double u3_u = log(2*x3_u + x1);

  long double u3 = u3_s*(u3_u-u3_l) + u3_l;

  long double x3 = 0.5*(exp(u3)- x1);


  ////////////////
  // Useful Definitions in Frame p1vec + qvec - p3vec = 0
  ////////////////
  long double Wc2 = 2*(x2-x3-x5) - x1 + SQR(ml2);

  long double E1 =  (x2 - x5)/sqrt(Wc2) ;
  long double E4 =  (Wc2+SQR(ml1))/2.0/sqrt(Wc2) ;
  long double q0 =  (x2-x1-x3)/sqrt(Wc2) ;
  long double qvec = sqrt(SQR(q0) + x1)  ;
  long double p4vec = (Wc2 - SQR(ml1))/2.0/sqrt(Wc2)  ;

  long double Cq = (q0*E1 - x2)/qvec/E1;
  long double Sq = sqrt(1 - SQR(Cq));

  long double E2 = (Wc2 - SQR(ml1))/2/sqrt(Wc2);


  ////////////////

  long double m6, u6, u6_l, u6_u, prop, jacob_u6;
  long double m6_l = 0;
  long double m6_u = 2*E1*E2;

  switch ((int) BSM) 
  {
    case (SMonly):
      u6_l = m6_l;
      u6_u = m6_u;
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6 = u6;
      jacob_u6 = 1.0;

      V2     = myMC->terms[3];
      A2     = myMC->terms[4];
      VA     = myMC->terms[5];
      break;

    case (INTERFERENCE):

      //////////////////////////////////////////////////////////
      // BSM PROPAGATOR 
      // long double prop = (1.0/(-2*(m6/mzprime/mzprime) - 1.0))/mzprime/mzprime;
      u6_l = log(1.0/(2.0*m6_u + mzprime*mzprime));
      u6_u = log(1.0/(2.0*m6_l + mzprime*mzprime));
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6   = -mzprime*mzprime/2.0 + exp(-u6)/2.0;

      prop = -exp(u6);
      jacob_u6 = exp(-u6)/2.0;

      V2     = myMC->terms[3] * prop;
      A2     = myMC->terms[4] * prop;
      VA     = myMC->terms[5] * prop;
      break;

    case (BSMonly):

      u6_l = 1.0/(2.0*m6_u + mzprime*mzprime);
      u6_u = 1.0/(2.0*m6_l + mzprime*mzprime);
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6   = -mzprime*mzprime/2.0 + 1.0/2.0/u6;

      prop = -u6;
      jacob_u6 = 1.0/(u6*u6)/2.0;

      V2     = myMC->terms[3] * prop*prop;
      A2     = myMC->terms[4] * prop*prop;
      VA     = myMC->terms[5] * prop*prop;
      break;
    
    case (SMandBSM):
     
      u6_l = m6_l;
      u6_u = m6_u;
      u6   = u6_s*(u6_u - u6_l) + u6_l;
      m6 = u6;
      jacob_u6 = 1.0;

      prop = 1.0/(2.0*m6 + mzprime*mzprime);

      V2     = (myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop) ;
      A2     = (myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop) ;
      VA     = (myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop);
      break;

    default:
      std::cout << "Could not determine what contributions (BSM) to include." << std::endl;
      break;
  }

  long double C2 = 1.0 - m6 / E1 / E2;
  long double S2 = sqrt(1.0 - SQR(C2));


  ////////////////
  long double PHI2 = PHI2_s*(2*M_PI-0.0) + 0.0;

  long double x4 = q0*E4 + (Sq*S2*cos(PHI2) + Cq*C2)*E2*qvec;

  long double C4 = (E4*q0 - x4)/(p4vec*qvec);
  long double S4 = sqrt(1.0 - SQR(C4));

  long double x6 = E1*E4 + E1*E2*C2; 

  //////////////////////////////////////////////////////////
  // Multiply by the appropriate Jacobians for our invariants and phase space factors
  long double Jacob =  (Wc2 - SQR(ml1))/8.0/Wc2/E1/E2/x2/pow(2*M_PI,6.0);

  // Jacobian due to change from CSW(+Matheus) to better integration variables AND VEGAS
  Jacob *=  (u1_u-u1_l)*exp(u1)*
            (u2_u-u2_l)*0.5*
            (u3_u-u3_l)*0.5*exp(u3)*
            (2*M_PI)*
            (x5_u-x5_l)*
            (u6_u-u6_l)*jacob_u6*
            (2*M_PI)*
            (2*M_PI);


  //////////////////////////////////////////////////////////
  // MATRIX ELEMENT ITSELF
  long double HH1, HH2;
  if (Z == 0)
  {
    HH1 = H1_n(sqrt(x1));
    HH2 = H2_n(sqrt(x1));
  }
  else if (Z == 1)
  {
    HH1 = H1_p(sqrt(x1));
    HH2 = H2_p(sqrt(x1));
  }
  else
  {
    std::cout<<"Error! Invalid Z. Diffractive regime has Z=0 or Z=1. "<<std::endl;
  }
  // includes Z*Z for the coherent enhancement
  long double dsigma = (4*(alphaQED*alphaQED)*(Gf*Gf)*(Diag22*((x1 + 2*x4)*(x1 + 2*x4))*(2*VA*(-2*(ml2*ml2)*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) - 2*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) + 2*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5)*x5 - 4*(x1*x1)*(x2*x2)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5 + 8*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*(x5*x5) - 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*(x5*x5*x5) - 3*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 4*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x6 + 4*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(-x1 + x2 - x3 - x4)*x6 + 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x5*x5)*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 - 2*(x1*x1)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6)*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x3*x5*(x2 - x5 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x3*x5*(x2 - x5 - x6)*x6 + 2*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) - A2*(2*(ml2*ml2)*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) + 2*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) - 2*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5)*x5 + 4*(x1*x1)*(x2*x2)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5 - 8*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*(x5*x5) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*(x5*x5*x5) + 3*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 4*ml1*(ml2*ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(x2 - x5 - x6) - 4*ml1*ml2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6) + 4*ml1*ml2*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6) + 4*ml1*ml2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*(x2 - x5 - x6) - 4*ml1*ml2*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x5*x5)*(x2 - x5 - x6) + 4*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x6 + 4*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(-x1 + x2 - x3 - x4)*x6 + 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x5*x5)*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 - 2*(x1*x1)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6)*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x3*x5*(x2 - x5 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x3*x5*(x2 - x5 - x6)*x6 + 2*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) - V2*(2*(ml2*ml2)*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) + 2*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5) - 2*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-(ml1*ml1) + ml2*ml2 - x1 + 2*x2 - 2*x3 - 2*x5)*x5 + 4*(x1*x1)*(x2*x2)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5 - 8*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*(x5*x5) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*(x5*x5*x5) + 3*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5*x5)*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - 4*ml1*(ml2*ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(x2 - x5 - x6) + 4*ml1*ml2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6) - 4*ml1*ml2*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6) - 4*ml1*ml2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*(x2 - x5 - x6) + 4*ml1*ml2*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x5*x5)*(x2 - x5 - x6) + 4*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x6 + 4*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x3*(-x1 + x2 - x3 - x4)*x6 + 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*(x5*x5)*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x5*x5)*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 - 2*(x1*x1)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6)*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x3*x5*(x2 - x5 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x3*x5*(x2 - x5 - x6)*x6 + 2*(ml2*ml2)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x5*x5)*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6))) + Diag11*((x1 + 2*x3)*(x1 + 2*x3))*(2*VA*(-4*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x5 - 4*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x4*x5 - 2*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) - x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*(ml1*ml1)*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) + 2*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x4*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) - 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 + 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*x6 + 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5*x6 + 2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 - 2*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6)*x6 + 4*(x1*x1)*(x2*x2)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 + 2*(x1*x1)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6)*x6 - 2*(x1*x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6)*x6 + 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x4*x5*(x2 - x5 - x6)*x6 - 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x4*x5*(x2 - x5 - x6)*x6 + 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) - 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*(x6*x6) - 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5*(x6*x6) - 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) - 8*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*(x6*x6) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*(x6*x6*x6) + 3*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) - A2*(4*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x5 + 4*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x4*x5 + 2*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*(ml1*ml1)*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) + 2*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x4*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) + 4*(ml1*ml1*ml1)*ml2*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x4*(x2 - x5 - x6) + 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5*x6 - 2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 - 2*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6)*x6 + 4*(x1*x1)*(x2*x2)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 - 4*ml1*ml2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6)*x6 + 4*ml1*ml2*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x2 - x5 - x6)*x6 - 2*(x1*x1)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6)*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x4*x5*(x2 - x5 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x4*x5*(x2 - x5 - x6)*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*(x6*x6) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5*(x6*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) - 8*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*(x6*x6) + 4*ml1*ml2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6)*(x6*x6) - 4*ml1*ml2*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x2 - x5 - x6)*(x6*x6) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*(x6*x6*x6) + 3*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6)) - V2*(4*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x5 + 4*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(-x1 + x2 - x3 - x4)*x4*x5 + 2*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5) + 2*(ml1*ml1)*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) + 2*(x2*x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x4*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6) - 4*(ml1*ml1*ml1)*ml2*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(x2 - x5 - x6) - 4*ml1*ml2*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6) + 4*ml1*ml2*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*x4*(x2 - x5 - x6) + 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*x6 - 4*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5*x6 - 2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*x6 - 2*(ml1*ml1)*(x2*x2)*(-2*(Enu*Enu)*HH1*(Mn*Mn)*x1 + 2*Enu*HH1*Mn*x1*x2 + (2*HH1*(Mn*Mn) - HH2*x1)*(x2*x2))*(ml1*ml1 - ml2*ml2 - x1 + 2*x2 - 2*x4 - 2*x6)*x6 + 4*(x1*x1)*(x2*x2)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*x6 + 4*ml1*ml2*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6)*x6 - 4*ml1*ml2*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x2 - x5 - x6)*x6 - 2*(x1*x1)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(x2 - x5 - x6)*x6 + 2*(x1*x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x5*(x2 - x5 - x6)*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x4*x5*(x2 - x5 - x6)*x6 + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*x4*x5*(x2 - x5 - x6)*x6 - 4*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x1 - x2 + x3 + x4)*x5*(x6*x6) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x1 - x2 + x3 + x4)*x5*(x6*x6) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((-(ml1*ml1) + ml2*ml2 - x1)/2. + x2 - x3 - x5)*x5*(x6*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x5*(ml1*ml1 - ml2*ml2 + x1 - 2*x2 + 2*x3 + 2*x5)*(x6*x6) - 8*(x1*x1)*x2*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*(x6*x6) - 4*ml1*ml2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x2 - x5 - x6)*(x6*x6) + 4*ml1*ml2*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*(x2 - x5 - x6)*(x6*x6) + 4*(x1*x1)*(HH2*(x2*x2) - HH1*((-2*Enu*Mn + x2)*(-2*Enu*Mn + x2)))*((ml1*ml1 - ml2*ml2 - x1)/2. + x2 - x4 - x6)*(x6*x6*x6) + 3*x1*(x2*x2)*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*x6*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) - 4*x1*x2*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6) + 2*x1*(2*(Enu*Enu)*HH1*(Mn*Mn)*x1 - 2*Enu*HH1*Mn*x1*x2 + (-2*HH1*(Mn*Mn) + HH2*x1)*(x2*x2))*(x6*x6*x6)*(-(ml1*ml1) + ml2*ml2 + x1 - 2*x2 + 2*x4 + 2*x6))) + 2*Diag12*(x1 + 2*x3)*(x1 + 2*x4)*(HH2*x1*(x2*x2*x2*x2)*(-2*V2*x1*(x2*x2) + 4*V2*(x2*x2*x2) + V2*x1*x2*x3 + 2*VA*x1*x2*x3 - 6*V2*(x2*x2)*x3 - 4*VA*(x2*x2)*x3 + 2*V2*x2*(x3*x3) + 4*VA*x2*(x3*x3) + V2*x1*x2*x4 - 2*VA*x1*x2*x4 - 6*V2*(x2*x2)*x4 + 4*VA*(x2*x2)*x4 + 4*V2*x2*x3*x4 + 2*V2*x2*(x4*x4) - 4*VA*x2*(x4*x4) + V2*(x1*x1)*x5 + 2*VA*(x1*x1)*x5 - 8*V2*(x2*x2)*x5 - 8*VA*(x2*x2)*x5 + 2*V2*x1*x3*x5 + 4*VA*x1*x3*x5 + 6*V2*x2*x3*x5 + 12*VA*x2*x3*x5 + 2*V2*x1*x4*x5 + 4*VA*x1*x4*x5 + 6*V2*x2*x4*x5 - 4*VA*x2*x4*x5 + 2*V2*x3*x4*x5 + 4*VA*x3*x4*x5 - 2*V2*(x4*x4)*x5 + 4*VA*(x4*x4)*x5 + 8*V2*x2*(x5*x5) + 16*VA*x2*(x5*x5) - 4*V2*x3*(x5*x5) - 8*VA*x3*(x5*x5) - 4*V2*(x5*x5*x5) - 8*VA*(x5*x5*x5) + V2*(x1*x1)*x6 - 2*VA*(x1*x1)*x6 - 8*V2*(x2*x2)*x6 + 8*VA*(x2*x2)*x6 + 2*V2*x1*x3*x6 - 4*VA*x1*x3*x6 + 6*V2*x2*x3*x6 + 4*VA*x2*x3*x6 - 2*V2*(x3*x3)*x6 - 4*VA*(x3*x3)*x6 + 2*V2*x1*x4*x6 - 4*VA*x1*x4*x6 + 6*V2*x2*x4*x6 - 12*VA*x2*x4*x6 + 2*V2*x3*x4*x6 - 4*VA*x3*x4*x6 + 4*V2*x1*x5*x6 + 8*V2*x2*x5*x6 - 4*V2*x3*x5*x6 - 8*VA*x3*x5*x6 - 4*V2*x4*x5*x6 + 8*VA*x4*x5*x6 - 4*V2*(x5*x5)*x6 - 8*VA*(x5*x5)*x6 + 8*V2*x2*(x6*x6) - 16*VA*x2*(x6*x6) - 4*V2*x4*(x6*x6) + 8*VA*x4*(x6*x6) - 4*V2*x5*(x6*x6) + 8*VA*x5*(x6*x6) - 4*V2*(x6*x6*x6) + 8*VA*(x6*x6*x6) + 2*(ml1*ml1*ml1)*ml2*V2*(-x2 + x5 + x6) - 2*ml1*ml2*V2*(-2*x1*x2 + 4*(x2*x2) - 3*x2*x3 - 3*x2*x4 - 4*x2*x5 + x3*x5 + x4*x5 + 2*(x5*x5) + ml2*ml2*(x2 - x5 - x6) - 4*x2*x6 + x3*x6 + x4*x6 + 4*x5*x6 + 2*(x6*x6)) + ml2*ml2*ml2*ml2*(2*VA*(x2 - x5 - x6) + V2*(-x5 + x6)) + ml1*ml1*ml1*ml1*(V2*(x5 - x6) + 2*VA*(-x2 + x5 + x6)) + ml1*ml1*(V2*(2*(x2*x2) - 2*x1*x5 + 2*x6*(-x3 + x4 + 2*x6) - x2*(x3 + x4 + 4*x6)) + 2*VA*(x1*x2 - 2*(x2*x2) - 2*(x5 + x6)*(x4 + 2*x6) + x2*(x3 + 3*x4 + 2*x5 + 6*x6))) - ml2*ml2*(2*VA*(x1*x2 - 2*(x2*x2) - 2*(x3 + 2*x5)*(x5 + x6) + x2*(3*x3 + x4 + 6*x5 + 2*x6)) + V2*(-2*(x2*x2) + x2*(x3 + x4 + 4*x5) + 2*(-(x3*x5) + x4*x5 - 2*(x5*x5) + x1*x6)))) + 2*(Enu*Enu)*HH1*(Mn*Mn)*x1*(-2*V2*x1*(x2*x2*x2*x2) + 4*V2*(x2*x2*x2*x2*x2) + V2*x1*(x2*x2*x2)*x3 + 2*VA*x1*(x2*x2*x2)*x3 - 6*V2*(x2*x2*x2*x2)*x3 - 4*VA*(x2*x2*x2*x2)*x3 + 2*V2*(x2*x2*x2)*(x3*x3) + 4*VA*(x2*x2*x2)*(x3*x3) + V2*x1*(x2*x2*x2)*x4 - 2*VA*x1*(x2*x2*x2)*x4 - 6*V2*(x2*x2*x2*x2)*x4 + 4*VA*(x2*x2*x2*x2)*x4 + 4*V2*(x2*x2*x2)*x3*x4 + 2*V2*(x2*x2*x2)*(x4*x4) - 4*VA*(x2*x2*x2)*(x4*x4) + 4*V2*(x1*x1)*(x2*x2)*x5 + 8*VA*(x1*x1)*(x2*x2)*x5 - 6*V2*x1*(x2*x2*x2)*x5 - 12*VA*x1*(x2*x2*x2)*x5 - 8*V2*(x2*x2*x2*x2)*x5 - 8*VA*(x2*x2*x2*x2)*x5 + 8*V2*x1*(x2*x2)*x3*x5 + 16*VA*x1*(x2*x2)*x3*x5 + 6*V2*(x2*x2*x2)*x3*x5 + 12*VA*(x2*x2*x2)*x3*x5 + 8*V2*x1*(x2*x2)*x4*x5 + 16*VA*x1*(x2*x2)*x4*x5 + 6*V2*(x2*x2*x2)*x4*x5 - 4*VA*(x2*x2*x2)*x4*x5 + 2*V2*(x2*x2)*x3*x4*x5 + 4*VA*(x2*x2)*x3*x4*x5 - 2*V2*(x2*x2)*(x4*x4)*x5 + 4*VA*(x2*x2)*(x4*x4)*x5 - 3*V2*(x1*x1)*x2*(x5*x5) - 6*VA*(x1*x1)*x2*(x5*x5) + 12*V2*x1*(x2*x2)*(x5*x5) + 24*VA*x1*(x2*x2)*(x5*x5) + 8*V2*(x2*x2*x2)*(x5*x5) + 16*VA*(x2*x2*x2)*(x5*x5) - 6*V2*x1*x2*x3*(x5*x5) - 12*VA*x1*x2*x3*(x5*x5) - 4*V2*(x2*x2)*x3*(x5*x5) - 8*VA*(x2*x2)*x3*(x5*x5) - 12*V2*x1*x2*x4*(x5*x5) - 24*VA*x1*x2*x4*(x5*x5) - 6*V2*x1*x2*(x5*x5*x5) - 12*VA*x1*x2*(x5*x5*x5) - 4*V2*(x2*x2)*(x5*x5*x5) - 8*VA*(x2*x2)*(x5*x5*x5) + 6*V2*x1*x4*(x5*x5*x5) + 12*VA*x1*x4*(x5*x5*x5) + 4*V2*(x1*x1)*(x2*x2)*x6 - 8*VA*(x1*x1)*(x2*x2)*x6 - 6*V2*x1*(x2*x2*x2)*x6 + 12*VA*x1*(x2*x2*x2)*x6 - 8*V2*(x2*x2*x2*x2)*x6 + 8*VA*(x2*x2*x2*x2)*x6 + 8*V2*x1*(x2*x2)*x3*x6 - 16*VA*x1*(x2*x2)*x3*x6 + 6*V2*(x2*x2*x2)*x3*x6 + 4*VA*(x2*x2*x2)*x3*x6 - 2*V2*(x2*x2)*(x3*x3)*x6 - 4*VA*(x2*x2)*(x3*x3)*x6 + 8*V2*x1*(x2*x2)*x4*x6 - 16*VA*x1*(x2*x2)*x4*x6 + 6*V2*(x2*x2*x2)*x4*x6 - 12*VA*(x2*x2*x2)*x4*x6 + 2*V2*(x2*x2)*x3*x4*x6 - 4*VA*(x2*x2)*x3*x4*x6 + 16*V2*x1*(x2*x2)*x5*x6 + 8*V2*(x2*x2*x2)*x5*x6 - 6*V2*x1*x2*x3*x5*x6 - 12*VA*x1*x2*x3*x5*x6 - 4*V2*(x2*x2)*x3*x5*x6 - 8*VA*(x2*x2)*x3*x5*x6 - 6*V2*x1*x2*x4*x5*x6 + 12*VA*x1*x2*x4*x5*x6 - 4*V2*(x2*x2)*x4*x5*x6 + 8*VA*(x2*x2)*x4*x5*x6 - 18*V2*x1*x2*(x5*x5)*x6 - 36*VA*x1*x2*(x5*x5)*x6 - 4*V2*(x2*x2)*(x5*x5)*x6 - 8*VA*(x2*x2)*(x5*x5)*x6 + 6*V2*x1*x3*(x5*x5)*x6 + 12*VA*x1*x3*(x5*x5)*x6 + 12*V2*x1*(x5*x5*x5)*x6 + 24*VA*x1*(x5*x5*x5)*x6 - 3*V2*(x1*x1)*x2*(x6*x6) + 6*VA*(x1*x1)*x2*(x6*x6) + 12*V2*x1*(x2*x2)*(x6*x6) - 24*VA*x1*(x2*x2)*(x6*x6) + 8*V2*(x2*x2*x2)*(x6*x6) - 16*VA*(x2*x2*x2)*(x6*x6) - 12*V2*x1*x2*x3*(x6*x6) + 24*VA*x1*x2*x3*(x6*x6) - 6*V2*x1*x2*x4*(x6*x6) + 12*VA*x1*x2*x4*(x6*x6) - 4*V2*(x2*x2)*x4*(x6*x6) + 8*VA*(x2*x2)*x4*(x6*x6) - 18*V2*x1*x2*x5*(x6*x6) + 36*VA*x1*x2*x5*(x6*x6) - 4*V2*(x2*x2)*x5*(x6*x6) + 8*VA*(x2*x2)*x5*(x6*x6) + 6*V2*x1*x4*x5*(x6*x6) - 12*VA*x1*x4*x5*(x6*x6) - 6*V2*x1*x2*(x6*x6*x6) + 12*VA*x1*x2*(x6*x6*x6) - 4*V2*(x2*x2)*(x6*x6*x6) + 8*VA*(x2*x2)*(x6*x6*x6) + 6*V2*x1*x3*(x6*x6*x6) - 12*VA*x1*x3*(x6*x6*x6) + 12*V2*x1*x5*(x6*x6*x6) - 24*VA*x1*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1)*ml2*V2*(x2*x2)*(-x2 + x5 + x6) + ml2*ml2*ml2*ml2*(x2*x2)*(2*VA*(x2 - x5 - x6) + V2*(-x5 + x6)) + ml1*ml1*ml1*ml1*(x2*x2)*(V2*(x5 - x6) + 2*VA*(-x2 + x5 + x6)) - 2*ml1*ml2*V2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6)))) - ml2*ml2*(V2*(-2*(x2*x2*x2*x2) + x2*x2*x2*(x3 + x4 + 4*x5) + x2*x2*(-2*x5*(x3 - x4 + 2*x5) + x1*(3*x5 - x6)) + 6*x1*x5*(x5 - x6)*x6 + 3*x1*x2*(-(x5*x5) + x6*x6)) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x3 + 2*x5)*(x5 + x6) + x2*(3*x3 + x4 + 6*x5 + 2*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)))) + ml1*ml1*(V2*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x5 + x6)*(x4 + 2*x6) + x2*(x3 + 3*x4 + 2*x5 + 6*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6))))) - 2*Enu*HH1*Mn*x1*x2*(-2*V2*x1*(x2*x2*x2*x2) + 4*V2*(x2*x2*x2*x2*x2) + V2*x1*(x2*x2*x2)*x3 + 2*VA*x1*(x2*x2*x2)*x3 - 6*V2*(x2*x2*x2*x2)*x3 - 4*VA*(x2*x2*x2*x2)*x3 + 2*V2*(x2*x2*x2)*(x3*x3) + 4*VA*(x2*x2*x2)*(x3*x3) + V2*x1*(x2*x2*x2)*x4 - 2*VA*x1*(x2*x2*x2)*x4 - 6*V2*(x2*x2*x2*x2)*x4 + 4*VA*(x2*x2*x2*x2)*x4 + 4*V2*(x2*x2*x2)*x3*x4 + 2*V2*(x2*x2*x2)*(x4*x4) - 4*VA*(x2*x2*x2)*(x4*x4) + 4*V2*(x1*x1)*(x2*x2)*x5 + 8*VA*(x1*x1)*(x2*x2)*x5 - 6*V2*x1*(x2*x2*x2)*x5 - 12*VA*x1*(x2*x2*x2)*x5 - 8*V2*(x2*x2*x2*x2)*x5 - 8*VA*(x2*x2*x2*x2)*x5 + 8*V2*x1*(x2*x2)*x3*x5 + 16*VA*x1*(x2*x2)*x3*x5 + 6*V2*(x2*x2*x2)*x3*x5 + 12*VA*(x2*x2*x2)*x3*x5 + 8*V2*x1*(x2*x2)*x4*x5 + 16*VA*x1*(x2*x2)*x4*x5 + 6*V2*(x2*x2*x2)*x4*x5 - 4*VA*(x2*x2*x2)*x4*x5 + 2*V2*(x2*x2)*x3*x4*x5 + 4*VA*(x2*x2)*x3*x4*x5 - 2*V2*(x2*x2)*(x4*x4)*x5 + 4*VA*(x2*x2)*(x4*x4)*x5 - 3*V2*(x1*x1)*x2*(x5*x5) - 6*VA*(x1*x1)*x2*(x5*x5) + 12*V2*x1*(x2*x2)*(x5*x5) + 24*VA*x1*(x2*x2)*(x5*x5) + 8*V2*(x2*x2*x2)*(x5*x5) + 16*VA*(x2*x2*x2)*(x5*x5) - 6*V2*x1*x2*x3*(x5*x5) - 12*VA*x1*x2*x3*(x5*x5) - 4*V2*(x2*x2)*x3*(x5*x5) - 8*VA*(x2*x2)*x3*(x5*x5) - 12*V2*x1*x2*x4*(x5*x5) - 24*VA*x1*x2*x4*(x5*x5) - 6*V2*x1*x2*(x5*x5*x5) - 12*VA*x1*x2*(x5*x5*x5) - 4*V2*(x2*x2)*(x5*x5*x5) - 8*VA*(x2*x2)*(x5*x5*x5) + 6*V2*x1*x4*(x5*x5*x5) + 12*VA*x1*x4*(x5*x5*x5) + 4*V2*(x1*x1)*(x2*x2)*x6 - 8*VA*(x1*x1)*(x2*x2)*x6 - 6*V2*x1*(x2*x2*x2)*x6 + 12*VA*x1*(x2*x2*x2)*x6 - 8*V2*(x2*x2*x2*x2)*x6 + 8*VA*(x2*x2*x2*x2)*x6 + 8*V2*x1*(x2*x2)*x3*x6 - 16*VA*x1*(x2*x2)*x3*x6 + 6*V2*(x2*x2*x2)*x3*x6 + 4*VA*(x2*x2*x2)*x3*x6 - 2*V2*(x2*x2)*(x3*x3)*x6 - 4*VA*(x2*x2)*(x3*x3)*x6 + 8*V2*x1*(x2*x2)*x4*x6 - 16*VA*x1*(x2*x2)*x4*x6 + 6*V2*(x2*x2*x2)*x4*x6 - 12*VA*(x2*x2*x2)*x4*x6 + 2*V2*(x2*x2)*x3*x4*x6 - 4*VA*(x2*x2)*x3*x4*x6 + 16*V2*x1*(x2*x2)*x5*x6 + 8*V2*(x2*x2*x2)*x5*x6 - 6*V2*x1*x2*x3*x5*x6 - 12*VA*x1*x2*x3*x5*x6 - 4*V2*(x2*x2)*x3*x5*x6 - 8*VA*(x2*x2)*x3*x5*x6 - 6*V2*x1*x2*x4*x5*x6 + 12*VA*x1*x2*x4*x5*x6 - 4*V2*(x2*x2)*x4*x5*x6 + 8*VA*(x2*x2)*x4*x5*x6 - 18*V2*x1*x2*(x5*x5)*x6 - 36*VA*x1*x2*(x5*x5)*x6 - 4*V2*(x2*x2)*(x5*x5)*x6 - 8*VA*(x2*x2)*(x5*x5)*x6 + 6*V2*x1*x3*(x5*x5)*x6 + 12*VA*x1*x3*(x5*x5)*x6 + 12*V2*x1*(x5*x5*x5)*x6 + 24*VA*x1*(x5*x5*x5)*x6 - 3*V2*(x1*x1)*x2*(x6*x6) + 6*VA*(x1*x1)*x2*(x6*x6) + 12*V2*x1*(x2*x2)*(x6*x6) - 24*VA*x1*(x2*x2)*(x6*x6) + 8*V2*(x2*x2*x2)*(x6*x6) - 16*VA*(x2*x2*x2)*(x6*x6) - 12*V2*x1*x2*x3*(x6*x6) + 24*VA*x1*x2*x3*(x6*x6) - 6*V2*x1*x2*x4*(x6*x6) + 12*VA*x1*x2*x4*(x6*x6) - 4*V2*(x2*x2)*x4*(x6*x6) + 8*VA*(x2*x2)*x4*(x6*x6) - 18*V2*x1*x2*x5*(x6*x6) + 36*VA*x1*x2*x5*(x6*x6) - 4*V2*(x2*x2)*x5*(x6*x6) + 8*VA*(x2*x2)*x5*(x6*x6) + 6*V2*x1*x4*x5*(x6*x6) - 12*VA*x1*x4*x5*(x6*x6) - 6*V2*x1*x2*(x6*x6*x6) + 12*VA*x1*x2*(x6*x6*x6) - 4*V2*(x2*x2)*(x6*x6*x6) + 8*VA*(x2*x2)*(x6*x6*x6) + 6*V2*x1*x3*(x6*x6*x6) - 12*VA*x1*x3*(x6*x6*x6) + 12*V2*x1*x5*(x6*x6*x6) - 24*VA*x1*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1)*ml2*V2*(x2*x2)*(-x2 + x5 + x6) + ml2*ml2*ml2*ml2*(x2*x2)*(2*VA*(x2 - x5 - x6) + V2*(-x5 + x6)) + ml1*ml1*ml1*ml1*(x2*x2)*(V2*(x5 - x6) + 2*VA*(-x2 + x5 + x6)) - 2*ml1*ml2*V2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6)))) - ml2*ml2*(V2*(-2*(x2*x2*x2*x2) + x2*x2*x2*(x3 + x4 + 4*x5) + x2*x2*(-2*x5*(x3 - x4 + 2*x5) + x1*(3*x5 - x6)) + 6*x1*x5*(x5 - x6)*x6 + 3*x1*x2*(-(x5*x5) + x6*x6)) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x3 + 2*x5)*(x5 + x6) + x2*(3*x3 + x4 + 6*x5 + 2*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)))) + ml1*ml1*(V2*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x5 + x6)*(x4 + 2*x6) + x2*(x3 + 3*x4 + 2*x5 + 6*x6)) + x1*(x2*x2*x2 + 3*(x2*x2)*(x5 + x6) + 6*x5*x6*(x5 + x6) - 3*x2*(x5*x5 + 4*x5*x6 + x6*x6))))) + HH1*(x2*x2)*(4*(Mn*Mn)*V2*x1*(x2*x2*x2*x2) - 8*(Mn*Mn)*V2*(x2*x2*x2*x2*x2) - 2*(Mn*Mn)*V2*x1*(x2*x2*x2)*x3 - 4*(Mn*Mn)*VA*x1*(x2*x2*x2)*x3 + 12*(Mn*Mn)*V2*(x2*x2*x2*x2)*x3 + 8*(Mn*Mn)*VA*(x2*x2*x2*x2)*x3 - 4*(Mn*Mn)*V2*(x2*x2*x2)*(x3*x3) - 8*(Mn*Mn)*VA*(x2*x2*x2)*(x3*x3) - 2*(Mn*Mn)*V2*x1*(x2*x2*x2)*x4 + 4*(Mn*Mn)*VA*x1*(x2*x2*x2)*x4 + 12*(Mn*Mn)*V2*(x2*x2*x2*x2)*x4 - 8*(Mn*Mn)*VA*(x2*x2*x2*x2)*x4 - 8*(Mn*Mn)*V2*(x2*x2*x2)*x3*x4 - 4*(Mn*Mn)*V2*(x2*x2*x2)*(x4*x4) + 8*(Mn*Mn)*VA*(x2*x2*x2)*(x4*x4) - 4*(Mn*Mn)*V2*(x1*x1)*(x2*x2)*x5 - 8*(Mn*Mn)*VA*(x1*x1)*(x2*x2)*x5 + V2*(x1*x1*x1)*(x2*x2)*x5 + 2*VA*(x1*x1*x1)*(x2*x2)*x5 + 4*(Mn*Mn)*V2*x1*(x2*x2*x2)*x5 + 8*(Mn*Mn)*VA*x1*(x2*x2*x2)*x5 - 2*V2*(x1*x1)*(x2*x2*x2)*x5 - 4*VA*(x1*x1)*(x2*x2*x2)*x5 + 16*(Mn*Mn)*V2*(x2*x2*x2*x2)*x5 + 16*(Mn*Mn)*VA*(x2*x2*x2*x2)*x5 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x3*x5 - 16*(Mn*Mn)*VA*x1*(x2*x2)*x3*x5 + 2*V2*(x1*x1)*(x2*x2)*x3*x5 + 4*VA*(x1*x1)*(x2*x2)*x3*x5 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x3*x5 - 24*(Mn*Mn)*VA*(x2*x2*x2)*x3*x5 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x4*x5 - 16*(Mn*Mn)*VA*x1*(x2*x2)*x4*x5 + 2*V2*(x1*x1)*(x2*x2)*x4*x5 + 4*VA*(x1*x1)*(x2*x2)*x4*x5 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x4*x5 + 8*(Mn*Mn)*VA*(x2*x2*x2)*x4*x5 - 4*(Mn*Mn)*V2*(x2*x2)*x3*x4*x5 - 8*(Mn*Mn)*VA*(x2*x2)*x3*x4*x5 + 4*(Mn*Mn)*V2*(x2*x2)*(x4*x4)*x5 - 8*(Mn*Mn)*VA*(x2*x2)*(x4*x4)*x5 + 2*(Mn*Mn)*V2*(x1*x1)*x2*(x5*x5) + 4*(Mn*Mn)*VA*(x1*x1)*x2*(x5*x5) - V2*(x1*x1*x1)*x2*(x5*x5) - 2*VA*(x1*x1*x1)*x2*(x5*x5) - 8*(Mn*Mn)*V2*x1*(x2*x2)*(x5*x5) - 16*(Mn*Mn)*VA*x1*(x2*x2)*(x5*x5) + 4*V2*(x1*x1)*(x2*x2)*(x5*x5) + 8*VA*(x1*x1)*(x2*x2)*(x5*x5) - 16*(Mn*Mn)*V2*(x2*x2*x2)*(x5*x5) - 32*(Mn*Mn)*VA*(x2*x2*x2)*(x5*x5) + 4*(Mn*Mn)*V2*x1*x2*x3*(x5*x5) + 8*(Mn*Mn)*VA*x1*x2*x3*(x5*x5) - 2*V2*(x1*x1)*x2*x3*(x5*x5) - 4*VA*(x1*x1)*x2*x3*(x5*x5) + 8*(Mn*Mn)*V2*(x2*x2)*x3*(x5*x5) + 16*(Mn*Mn)*VA*(x2*x2)*x3*(x5*x5) + 8*(Mn*Mn)*V2*x1*x2*x4*(x5*x5) + 16*(Mn*Mn)*VA*x1*x2*x4*(x5*x5) - 4*V2*(x1*x1)*x2*x4*(x5*x5) - 8*VA*(x1*x1)*x2*x4*(x5*x5) + 4*(Mn*Mn)*V2*x1*x2*(x5*x5*x5) + 8*(Mn*Mn)*VA*x1*x2*(x5*x5*x5) - 2*V2*(x1*x1)*x2*(x5*x5*x5) - 4*VA*(x1*x1)*x2*(x5*x5*x5) + 8*(Mn*Mn)*V2*(x2*x2)*(x5*x5*x5) + 16*(Mn*Mn)*VA*(x2*x2)*(x5*x5*x5) - 4*(Mn*Mn)*V2*x1*x4*(x5*x5*x5) - 8*(Mn*Mn)*VA*x1*x4*(x5*x5*x5) + 2*V2*(x1*x1)*x4*(x5*x5*x5) + 4*VA*(x1*x1)*x4*(x5*x5*x5) + 4*(ml1*ml1*ml1)*ml2*(Mn*Mn)*V2*(x2*x2)*(x2 - x5 - x6) - 4*(Mn*Mn)*V2*(x1*x1)*(x2*x2)*x6 + 8*(Mn*Mn)*VA*(x1*x1)*(x2*x2)*x6 + V2*(x1*x1*x1)*(x2*x2)*x6 - 2*VA*(x1*x1*x1)*(x2*x2)*x6 + 4*(Mn*Mn)*V2*x1*(x2*x2*x2)*x6 - 8*(Mn*Mn)*VA*x1*(x2*x2*x2)*x6 - 2*V2*(x1*x1)*(x2*x2*x2)*x6 + 4*VA*(x1*x1)*(x2*x2*x2)*x6 + 16*(Mn*Mn)*V2*(x2*x2*x2*x2)*x6 - 16*(Mn*Mn)*VA*(x2*x2*x2*x2)*x6 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x3*x6 + 16*(Mn*Mn)*VA*x1*(x2*x2)*x3*x6 + 2*V2*(x1*x1)*(x2*x2)*x3*x6 - 4*VA*(x1*x1)*(x2*x2)*x3*x6 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x3*x6 - 8*(Mn*Mn)*VA*(x2*x2*x2)*x3*x6 + 4*(Mn*Mn)*V2*(x2*x2)*(x3*x3)*x6 + 8*(Mn*Mn)*VA*(x2*x2)*(x3*x3)*x6 - 8*(Mn*Mn)*V2*x1*(x2*x2)*x4*x6 + 16*(Mn*Mn)*VA*x1*(x2*x2)*x4*x6 + 2*V2*(x1*x1)*(x2*x2)*x4*x6 - 4*VA*(x1*x1)*(x2*x2)*x4*x6 - 12*(Mn*Mn)*V2*(x2*x2*x2)*x4*x6 + 24*(Mn*Mn)*VA*(x2*x2*x2)*x4*x6 - 4*(Mn*Mn)*V2*(x2*x2)*x3*x4*x6 + 8*(Mn*Mn)*VA*(x2*x2)*x3*x4*x6 - 16*(Mn*Mn)*V2*x1*(x2*x2)*x5*x6 + 4*V2*(x1*x1)*(x2*x2)*x5*x6 - 16*(Mn*Mn)*V2*(x2*x2*x2)*x5*x6 + 4*(Mn*Mn)*V2*x1*x2*x3*x5*x6 + 8*(Mn*Mn)*VA*x1*x2*x3*x5*x6 - 2*V2*(x1*x1)*x2*x3*x5*x6 - 4*VA*(x1*x1)*x2*x3*x5*x6 + 8*(Mn*Mn)*V2*(x2*x2)*x3*x5*x6 + 16*(Mn*Mn)*VA*(x2*x2)*x3*x5*x6 + 4*(Mn*Mn)*V2*x1*x2*x4*x5*x6 - 8*(Mn*Mn)*VA*x1*x2*x4*x5*x6 - 2*V2*(x1*x1)*x2*x4*x5*x6 + 4*VA*(x1*x1)*x2*x4*x5*x6 + 8*(Mn*Mn)*V2*(x2*x2)*x4*x5*x6 - 16*(Mn*Mn)*VA*(x2*x2)*x4*x5*x6 + 12*(Mn*Mn)*V2*x1*x2*(x5*x5)*x6 + 24*(Mn*Mn)*VA*x1*x2*(x5*x5)*x6 - 6*V2*(x1*x1)*x2*(x5*x5)*x6 - 12*VA*(x1*x1)*x2*(x5*x5)*x6 + 8*(Mn*Mn)*V2*(x2*x2)*(x5*x5)*x6 + 16*(Mn*Mn)*VA*(x2*x2)*(x5*x5)*x6 - 4*(Mn*Mn)*V2*x1*x3*(x5*x5)*x6 - 8*(Mn*Mn)*VA*x1*x3*(x5*x5)*x6 + 2*V2*(x1*x1)*x3*(x5*x5)*x6 + 4*VA*(x1*x1)*x3*(x5*x5)*x6 - 8*(Mn*Mn)*V2*x1*(x5*x5*x5)*x6 - 16*(Mn*Mn)*VA*x1*(x5*x5*x5)*x6 + 4*V2*(x1*x1)*(x5*x5*x5)*x6 + 8*VA*(x1*x1)*(x5*x5*x5)*x6 + 2*(Mn*Mn)*V2*(x1*x1)*x2*(x6*x6) - 4*(Mn*Mn)*VA*(x1*x1)*x2*(x6*x6) - V2*(x1*x1*x1)*x2*(x6*x6) + 2*VA*(x1*x1*x1)*x2*(x6*x6) - 8*(Mn*Mn)*V2*x1*(x2*x2)*(x6*x6) + 16*(Mn*Mn)*VA*x1*(x2*x2)*(x6*x6) + 4*V2*(x1*x1)*(x2*x2)*(x6*x6) - 8*VA*(x1*x1)*(x2*x2)*(x6*x6) - 16*(Mn*Mn)*V2*(x2*x2*x2)*(x6*x6) + 32*(Mn*Mn)*VA*(x2*x2*x2)*(x6*x6) + 8*(Mn*Mn)*V2*x1*x2*x3*(x6*x6) - 16*(Mn*Mn)*VA*x1*x2*x3*(x6*x6) - 4*V2*(x1*x1)*x2*x3*(x6*x6) + 8*VA*(x1*x1)*x2*x3*(x6*x6) + 4*(Mn*Mn)*V2*x1*x2*x4*(x6*x6) - 8*(Mn*Mn)*VA*x1*x2*x4*(x6*x6) - 2*V2*(x1*x1)*x2*x4*(x6*x6) + 4*VA*(x1*x1)*x2*x4*(x6*x6) + 8*(Mn*Mn)*V2*(x2*x2)*x4*(x6*x6) - 16*(Mn*Mn)*VA*(x2*x2)*x4*(x6*x6) + 12*(Mn*Mn)*V2*x1*x2*x5*(x6*x6) - 24*(Mn*Mn)*VA*x1*x2*x5*(x6*x6) - 6*V2*(x1*x1)*x2*x5*(x6*x6) + 12*VA*(x1*x1)*x2*x5*(x6*x6) + 8*(Mn*Mn)*V2*(x2*x2)*x5*(x6*x6) - 16*(Mn*Mn)*VA*(x2*x2)*x5*(x6*x6) - 4*(Mn*Mn)*V2*x1*x4*x5*(x6*x6) + 8*(Mn*Mn)*VA*x1*x4*x5*(x6*x6) + 2*V2*(x1*x1)*x4*x5*(x6*x6) - 4*VA*(x1*x1)*x4*x5*(x6*x6) + 4*(Mn*Mn)*V2*x1*x2*(x6*x6*x6) - 8*(Mn*Mn)*VA*x1*x2*(x6*x6*x6) - 2*V2*(x1*x1)*x2*(x6*x6*x6) + 4*VA*(x1*x1)*x2*(x6*x6*x6) + 8*(Mn*Mn)*V2*(x2*x2)*(x6*x6*x6) - 16*(Mn*Mn)*VA*(x2*x2)*(x6*x6*x6) - 4*(Mn*Mn)*V2*x1*x3*(x6*x6*x6) + 8*(Mn*Mn)*VA*x1*x3*(x6*x6*x6) + 2*V2*(x1*x1)*x3*(x6*x6*x6) - 4*VA*(x1*x1)*x3*(x6*x6*x6) - 8*(Mn*Mn)*V2*x1*x5*(x6*x6*x6) + 16*(Mn*Mn)*VA*x1*x5*(x6*x6*x6) + 4*V2*(x1*x1)*x5*(x6*x6*x6) - 8*VA*(x1*x1)*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1*ml1)*(Mn*Mn)*(x2*x2)*(2*VA*(x2 - x5 - x6) + V2*(-x5 + x6)) + 2*(ml2*ml2*ml2*ml2)*(Mn*Mn)*(x2*x2)*(V2*(x5 - x6) + 2*VA*(-x2 + x5 + x6)) + 2*ml1*ml2*V2*(2*(ml2*ml2)*(Mn*Mn)*(x2*x2)*(x2 - x5 - x6) - x1*x1*(x2*x2*x2 - 2*(x2*x2)*(x5 + x6) - 2*x5*x6*(x5 + x6) + x2*(x5*x5 + 4*x5*x6 + x6*x6)) + 2*(Mn*Mn)*(-(x1*(x2*x2*x2 + 2*(x2*x2)*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6))))) + ml2*ml2*(x1*x1*(-(V2*(x5 - x6)*(x2*x2 + 2*x5*x6 - x2*(x5 + x6))) - 2*VA*(x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) + 2*(Mn*Mn)*(V2*(-2*(x2*x2*x2*x2) + x2*x2*x2*(x3 + x4 + 4*x5) + 2*x1*x5*(x5 - x6)*x6 + x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(-2*x5*(x3 - x4 + 2*x5) + x1*(x5 + x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x3 + 2*x5)*(x5 + x6) + x2*(3*x3 + x4 + 6*x5 + 2*x6)) + x1*(x2*x2*x2 + x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))))) + ml1*ml1*(x1*x1*(V2*(x5 - x6)*(x2*x2 + 2*x5*x6 - x2*(x5 + x6)) + 2*VA*(x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) - 2*(Mn*Mn)*(V2*(2*(x2*x2*x2*x2) + 2*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + x1*x2*(-(x5*x5) + x6*x6) - x2*x2*(2*(x3 - x4 - 2*x6)*x6 + x1*(x5 + x6))) + 2*VA*(x2*x2*(-2*(x2*x2) - 2*(x5 + x6)*(x4 + 2*x6) + x2*(x3 + 3*x4 + 2*x5 + 6*x6)) + x1*(x2*x2*x2 + x2*x2*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6)))))) + A2*(HH2*x1*(x2*x2*x2*x2)*(-2*x1*(x2*x2) + 4*(x2*x2*x2) + x1*x2*x3 - 6*(x2*x2)*x3 + 2*x2*(x3*x3) + x1*x2*x4 - 6*(x2*x2)*x4 + 4*x2*x3*x4 + 2*x2*(x4*x4) + x1*x1*x5 - 8*(x2*x2)*x5 + 2*x1*x3*x5 + 6*x2*x3*x5 + 2*x1*x4*x5 + 6*x2*x4*x5 + 2*x3*x4*x5 - 2*(x4*x4)*x5 + 8*x2*(x5*x5) - 4*x3*(x5*x5) - 4*(x5*x5*x5) + 2*(ml1*ml1*ml1)*ml2*(x2 - x5 - x6) + ml1*ml1*ml1*ml1*(x5 - x6) + x1*x1*x6 - 8*(x2*x2)*x6 + 2*x1*x3*x6 + 6*x2*x3*x6 - 2*(x3*x3)*x6 + 2*x1*x4*x6 + 6*x2*x4*x6 + 2*x3*x4*x6 + 4*x1*x5*x6 + 8*x2*x5*x6 - 4*x3*x5*x6 - 4*x4*x5*x6 - 4*(x5*x5)*x6 + 8*x2*(x6*x6) - 4*x4*(x6*x6) - 4*x5*(x6*x6) - 4*(x6*x6*x6) + ml2*ml2*ml2*ml2*(-x5 + x6) + ml2*ml2*(2*(x2*x2) + 2*x3*x5 - 2*x4*x5 + 4*(x5*x5) - x2*(x3 + x4 + 4*x5) - 2*x1*x6) + 2*ml1*ml2*(-2*x1*x2 + 4*(x2*x2) - 3*x2*x3 - 3*x2*x4 - 4*x2*x5 + x3*x5 + x4*x5 + 2*(x5*x5) + ml2*ml2*(x2 - x5 - x6) - 4*x2*x6 + x3*x6 + x4*x6 + 4*x5*x6 + 2*(x6*x6)) + ml1*ml1*(2*(x2*x2) - 2*x1*x5 + 2*x6*(-x3 + x4 + 2*x6) - x2*(x3 + x4 + 4*x6))) + 2*(Enu*Enu)*HH1*(Mn*Mn)*x1*(-2*x1*(x2*x2*x2*x2) + 4*(x2*x2*x2*x2*x2) + x1*(x2*x2*x2)*x3 - 6*(x2*x2*x2*x2)*x3 + 2*(x2*x2*x2)*(x3*x3) + x1*(x2*x2*x2)*x4 - 6*(x2*x2*x2*x2)*x4 + 4*(x2*x2*x2)*x3*x4 + 2*(x2*x2*x2)*(x4*x4) + 4*(x1*x1)*(x2*x2)*x5 - 6*x1*(x2*x2*x2)*x5 - 8*(x2*x2*x2*x2)*x5 + 8*x1*(x2*x2)*x3*x5 + 6*(x2*x2*x2)*x3*x5 + 8*x1*(x2*x2)*x4*x5 + 6*(x2*x2*x2)*x4*x5 + 2*(x2*x2)*x3*x4*x5 - 2*(x2*x2)*(x4*x4)*x5 - 3*(x1*x1)*x2*(x5*x5) + 12*x1*(x2*x2)*(x5*x5) + 8*(x2*x2*x2)*(x5*x5) - 6*x1*x2*x3*(x5*x5) - 4*(x2*x2)*x3*(x5*x5) - 12*x1*x2*x4*(x5*x5) - 6*x1*x2*(x5*x5*x5) - 4*(x2*x2)*(x5*x5*x5) + 6*x1*x4*(x5*x5*x5) + 2*(ml1*ml1*ml1)*ml2*(x2*x2)*(x2 - x5 - x6) + ml1*ml1*ml1*ml1*(x2*x2)*(x5 - x6) + 4*(x1*x1)*(x2*x2)*x6 - 6*x1*(x2*x2*x2)*x6 - 8*(x2*x2*x2*x2)*x6 + 8*x1*(x2*x2)*x3*x6 + 6*(x2*x2*x2)*x3*x6 - 2*(x2*x2)*(x3*x3)*x6 + 8*x1*(x2*x2)*x4*x6 + 6*(x2*x2*x2)*x4*x6 + 2*(x2*x2)*x3*x4*x6 + 16*x1*(x2*x2)*x5*x6 + 8*(x2*x2*x2)*x5*x6 - 6*x1*x2*x3*x5*x6 - 4*(x2*x2)*x3*x5*x6 - 6*x1*x2*x4*x5*x6 - 4*(x2*x2)*x4*x5*x6 - 18*x1*x2*(x5*x5)*x6 - 4*(x2*x2)*(x5*x5)*x6 + 6*x1*x3*(x5*x5)*x6 + 12*x1*(x5*x5*x5)*x6 - 3*(x1*x1)*x2*(x6*x6) + 12*x1*(x2*x2)*(x6*x6) + 8*(x2*x2*x2)*(x6*x6) - 12*x1*x2*x3*(x6*x6) - 6*x1*x2*x4*(x6*x6) - 4*(x2*x2)*x4*(x6*x6) - 18*x1*x2*x5*(x6*x6) - 4*(x2*x2)*x5*(x6*x6) + 6*x1*x4*x5*(x6*x6) - 6*x1*x2*(x6*x6*x6) - 4*(x2*x2)*(x6*x6*x6) + 6*x1*x3*(x6*x6*x6) + 12*x1*x5*(x6*x6*x6) + ml2*ml2*ml2*ml2*(x2*x2)*(-x5 + x6) + ml2*ml2*(2*(x2*x2*x2*x2) - x2*x2*x2*(x3 + x4 + 4*x5) + 6*x1*x5*x6*(-x5 + x6) + 3*x1*x2*(x5*x5 - x6*x6) + x2*x2*(2*x5*(x3 - x4 + 2*x5) + x1*(-3*x5 + x6))) + ml1*ml1*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*ml1*ml2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6))))) - 2*Enu*HH1*Mn*x1*x2*(-2*x1*(x2*x2*x2*x2) + 4*(x2*x2*x2*x2*x2) + x1*(x2*x2*x2)*x3 - 6*(x2*x2*x2*x2)*x3 + 2*(x2*x2*x2)*(x3*x3) + x1*(x2*x2*x2)*x4 - 6*(x2*x2*x2*x2)*x4 + 4*(x2*x2*x2)*x3*x4 + 2*(x2*x2*x2)*(x4*x4) + 4*(x1*x1)*(x2*x2)*x5 - 6*x1*(x2*x2*x2)*x5 - 8*(x2*x2*x2*x2)*x5 + 8*x1*(x2*x2)*x3*x5 + 6*(x2*x2*x2)*x3*x5 + 8*x1*(x2*x2)*x4*x5 + 6*(x2*x2*x2)*x4*x5 + 2*(x2*x2)*x3*x4*x5 - 2*(x2*x2)*(x4*x4)*x5 - 3*(x1*x1)*x2*(x5*x5) + 12*x1*(x2*x2)*(x5*x5) + 8*(x2*x2*x2)*(x5*x5) - 6*x1*x2*x3*(x5*x5) - 4*(x2*x2)*x3*(x5*x5) - 12*x1*x2*x4*(x5*x5) - 6*x1*x2*(x5*x5*x5) - 4*(x2*x2)*(x5*x5*x5) + 6*x1*x4*(x5*x5*x5) + 2*(ml1*ml1*ml1)*ml2*(x2*x2)*(x2 - x5 - x6) + ml1*ml1*ml1*ml1*(x2*x2)*(x5 - x6) + 4*(x1*x1)*(x2*x2)*x6 - 6*x1*(x2*x2*x2)*x6 - 8*(x2*x2*x2*x2)*x6 + 8*x1*(x2*x2)*x3*x6 + 6*(x2*x2*x2)*x3*x6 - 2*(x2*x2)*(x3*x3)*x6 + 8*x1*(x2*x2)*x4*x6 + 6*(x2*x2*x2)*x4*x6 + 2*(x2*x2)*x3*x4*x6 + 16*x1*(x2*x2)*x5*x6 + 8*(x2*x2*x2)*x5*x6 - 6*x1*x2*x3*x5*x6 - 4*(x2*x2)*x3*x5*x6 - 6*x1*x2*x4*x5*x6 - 4*(x2*x2)*x4*x5*x6 - 18*x1*x2*(x5*x5)*x6 - 4*(x2*x2)*(x5*x5)*x6 + 6*x1*x3*(x5*x5)*x6 + 12*x1*(x5*x5*x5)*x6 - 3*(x1*x1)*x2*(x6*x6) + 12*x1*(x2*x2)*(x6*x6) + 8*(x2*x2*x2)*(x6*x6) - 12*x1*x2*x3*(x6*x6) - 6*x1*x2*x4*(x6*x6) - 4*(x2*x2)*x4*(x6*x6) - 18*x1*x2*x5*(x6*x6) - 4*(x2*x2)*x5*(x6*x6) + 6*x1*x4*x5*(x6*x6) - 6*x1*x2*(x6*x6*x6) - 4*(x2*x2)*(x6*x6*x6) + 6*x1*x3*(x6*x6*x6) + 12*x1*x5*(x6*x6*x6) + ml2*ml2*ml2*ml2*(x2*x2)*(-x5 + x6) + ml2*ml2*(2*(x2*x2*x2*x2) - x2*x2*x2*(x3 + x4 + 4*x5) + 6*x1*x5*x6*(-x5 + x6) + 3*x1*x2*(x5*x5 - x6*x6) + x2*x2*(2*x5*(x3 - x4 + 2*x5) + x1*(-3*x5 + x6))) + ml1*ml1*(2*(x2*x2*x2*x2) + 6*x1*x5*(x5 - x6)*x6 - x2*x2*x2*(x3 + x4 + 4*x6) + 3*x1*x2*(-(x5*x5) + x6*x6) + x2*x2*(x1*(x5 - 3*x6) + 2*x6*(-x3 + x4 + 2*x6))) + 2*ml1*ml2*(ml2*ml2*(x2*x2)*(x2 - x5 - x6) + x1*(x2*x2*x2 - 6*(x2*x2)*(x5 + x6) - 6*x5*x6*(x5 + x6) + 3*x2*(x5*x5 + 4*x5*x6 + x6*x6)) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6))))) + HH1*(x2*x2)*(-4*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2*x2) + 4*(Mn*Mn)*x1*(x2*x2*x2*x2) - 8*(Mn*Mn)*(x2*x2*x2*x2*x2) + 2*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2)*x3 - 2*(Mn*Mn)*x1*(x2*x2*x2)*x3 + 12*(Mn*Mn)*(x2*x2*x2*x2)*x3 - 4*(Mn*Mn)*(x2*x2*x2)*(x3*x3) + 2*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2)*x4 - 2*(Mn*Mn)*x1*(x2*x2*x2)*x4 + 12*(Mn*Mn)*(x2*x2*x2*x2)*x4 - 8*(Mn*Mn)*(x2*x2*x2)*x3*x4 - 4*(Mn*Mn)*(x2*x2*x2)*(x4*x4) + 2*(ml2*ml2*ml2*ml2)*(Mn*Mn)*(x2*x2)*x5 + 2*(ml2*ml2)*(Mn*Mn)*x1*(x2*x2)*x5 - ml2*ml2*(x1*x1)*(x2*x2)*x5 - 4*(Mn*Mn)*(x1*x1)*(x2*x2)*x5 + x1*x1*x1*(x2*x2)*x5 + 8*(ml2*ml2)*(Mn*Mn)*(x2*x2*x2)*x5 + 4*(Mn*Mn)*x1*(x2*x2*x2)*x5 - 2*(x1*x1)*(x2*x2*x2)*x5 + 16*(Mn*Mn)*(x2*x2*x2*x2)*x5 - 4*(ml2*ml2)*(Mn*Mn)*(x2*x2)*x3*x5 - 8*(Mn*Mn)*x1*(x2*x2)*x3*x5 + 2*(x1*x1)*(x2*x2)*x3*x5 - 12*(Mn*Mn)*(x2*x2*x2)*x3*x5 + 4*(ml2*ml2)*(Mn*Mn)*(x2*x2)*x4*x5 - 8*(Mn*Mn)*x1*(x2*x2)*x4*x5 + 2*(x1*x1)*(x2*x2)*x4*x5 - 12*(Mn*Mn)*(x2*x2*x2)*x4*x5 - 4*(Mn*Mn)*(x2*x2)*x3*x4*x5 + 4*(Mn*Mn)*(x2*x2)*(x4*x4)*x5 - 2*(ml2*ml2)*(Mn*Mn)*x1*x2*(x5*x5) + ml2*ml2*(x1*x1)*x2*(x5*x5) + 2*(Mn*Mn)*(x1*x1)*x2*(x5*x5) - x1*x1*x1*x2*(x5*x5) - 8*(ml2*ml2)*(Mn*Mn)*(x2*x2)*(x5*x5) - 8*(Mn*Mn)*x1*(x2*x2)*(x5*x5) + 4*(x1*x1)*(x2*x2)*(x5*x5) - 16*(Mn*Mn)*(x2*x2*x2)*(x5*x5) + 4*(Mn*Mn)*x1*x2*x3*(x5*x5) - 2*(x1*x1)*x2*x3*(x5*x5) + 8*(Mn*Mn)*(x2*x2)*x3*(x5*x5) + 8*(Mn*Mn)*x1*x2*x4*(x5*x5) - 4*(x1*x1)*x2*x4*(x5*x5) + 4*(Mn*Mn)*x1*x2*(x5*x5*x5) - 2*(x1*x1)*x2*(x5*x5*x5) + 8*(Mn*Mn)*(x2*x2)*(x5*x5*x5) - 4*(Mn*Mn)*x1*x4*(x5*x5*x5) + 2*(x1*x1)*x4*(x5*x5*x5) - 2*(ml2*ml2*ml2*ml2)*(Mn*Mn)*(x2*x2)*x6 + 2*(ml2*ml2)*(Mn*Mn)*x1*(x2*x2)*x6 + ml2*ml2*(x1*x1)*(x2*x2)*x6 - 4*(Mn*Mn)*(x1*x1)*(x2*x2)*x6 + x1*x1*x1*(x2*x2)*x6 + 4*(Mn*Mn)*x1*(x2*x2*x2)*x6 - 2*(x1*x1)*(x2*x2*x2)*x6 + 16*(Mn*Mn)*(x2*x2*x2*x2)*x6 - 8*(Mn*Mn)*x1*(x2*x2)*x3*x6 + 2*(x1*x1)*(x2*x2)*x3*x6 - 12*(Mn*Mn)*(x2*x2*x2)*x3*x6 + 4*(Mn*Mn)*(x2*x2)*(x3*x3)*x6 - 8*(Mn*Mn)*x1*(x2*x2)*x4*x6 + 2*(x1*x1)*(x2*x2)*x4*x6 - 12*(Mn*Mn)*(x2*x2*x2)*x4*x6 - 4*(Mn*Mn)*(x2*x2)*x3*x4*x6 - 16*(Mn*Mn)*x1*(x2*x2)*x5*x6 + 4*(x1*x1)*(x2*x2)*x5*x6 - 16*(Mn*Mn)*(x2*x2*x2)*x5*x6 + 4*(Mn*Mn)*x1*x2*x3*x5*x6 - 2*(x1*x1)*x2*x3*x5*x6 + 8*(Mn*Mn)*(x2*x2)*x3*x5*x6 + 4*(Mn*Mn)*x1*x2*x4*x5*x6 - 2*(x1*x1)*x2*x4*x5*x6 + 8*(Mn*Mn)*(x2*x2)*x4*x5*x6 + 4*(ml2*ml2)*(Mn*Mn)*x1*(x5*x5)*x6 - 2*(ml2*ml2)*(x1*x1)*(x5*x5)*x6 + 12*(Mn*Mn)*x1*x2*(x5*x5)*x6 - 6*(x1*x1)*x2*(x5*x5)*x6 + 8*(Mn*Mn)*(x2*x2)*(x5*x5)*x6 - 4*(Mn*Mn)*x1*x3*(x5*x5)*x6 + 2*(x1*x1)*x3*(x5*x5)*x6 - 8*(Mn*Mn)*x1*(x5*x5*x5)*x6 + 4*(x1*x1)*(x5*x5*x5)*x6 + 2*(ml2*ml2)*(Mn*Mn)*x1*x2*(x6*x6) - ml2*ml2*(x1*x1)*x2*(x6*x6) + 2*(Mn*Mn)*(x1*x1)*x2*(x6*x6) - x1*x1*x1*x2*(x6*x6) - 8*(Mn*Mn)*x1*(x2*x2)*(x6*x6) + 4*(x1*x1)*(x2*x2)*(x6*x6) - 16*(Mn*Mn)*(x2*x2*x2)*(x6*x6) + 8*(Mn*Mn)*x1*x2*x3*(x6*x6) - 4*(x1*x1)*x2*x3*(x6*x6) + 4*(Mn*Mn)*x1*x2*x4*(x6*x6) - 2*(x1*x1)*x2*x4*(x6*x6) + 8*(Mn*Mn)*(x2*x2)*x4*(x6*x6) - 4*(ml2*ml2)*(Mn*Mn)*x1*x5*(x6*x6) + 2*(ml2*ml2)*(x1*x1)*x5*(x6*x6) + 12*(Mn*Mn)*x1*x2*x5*(x6*x6) - 6*(x1*x1)*x2*x5*(x6*x6) + 8*(Mn*Mn)*(x2*x2)*x5*(x6*x6) - 4*(Mn*Mn)*x1*x4*x5*(x6*x6) + 2*(x1*x1)*x4*x5*(x6*x6) + 4*(Mn*Mn)*x1*x2*(x6*x6*x6) - 2*(x1*x1)*x2*(x6*x6*x6) + 8*(Mn*Mn)*(x2*x2)*(x6*x6*x6) - 4*(Mn*Mn)*x1*x3*(x6*x6*x6) + 2*(x1*x1)*x3*(x6*x6*x6) - 8*(Mn*Mn)*x1*x5*(x6*x6*x6) + 4*(x1*x1)*x5*(x6*x6*x6) + 2*(ml1*ml1*ml1*ml1)*(Mn*Mn)*(x2*x2)*(-x5 + x6) + 4*(ml1*ml1*ml1)*ml2*(Mn*Mn)*(x2*x2)*(-x2 + x5 + x6) + ml1*ml1*(x1*x1*(x5 - x6)*(x2*x2 + 2*x5*x6 - x2*(x5 + x6)) + 2*(Mn*Mn)*(-2*(x2*x2*x2*x2) + 2*x1*x5*x6*(-x5 + x6) + x2*x2*x2*(x3 + x4 + 4*x6) + x1*x2*(x5*x5 - x6*x6) + x2*x2*(2*(x3 - x4 - 2*x6)*x6 + x1*(x5 + x6)))) - 2*ml1*ml2*(2*(ml2*ml2)*(Mn*Mn)*(x2*x2)*(x2 - x5 - x6) - x1*x1*(x2*x2*x2 - 2*(x2*x2)*(x5 + x6) - 2*x5*x6*(x5 + x6) + x2*(x5*x5 + 4*x5*x6 + x6*x6)) + 2*(Mn*Mn)*(-(x1*(x2*x2*x2 + 2*(x2*x2)*(x5 + x6) + 2*x5*x6*(x5 + x6) - x2*(x5*x5 + 4*x5*x6 + x6*x6))) + x2*x2*(4*(x2*x2) + (x5 + x6)*(x3 + x4 + 2*(x5 + x6)) - x2*(3*x3 + 3*x4 + 4*(x5 + x6))))))))))/(Enu*Enu*(Mn*Mn)*(x1*x1)*(x2*x2*x2*x2)*((x1 + 2*x3)*(x1 + 2*x3))*((x1 + 2*x4)*(x1 + 2*x4)));

  // FERMI BLOCKING
  if(myMC->PAULI_BLOCKING == W_BLOCKING){

    dsigma *= F_PAULI_BLOCK(x1);
  }

  ///////////////////////////////////////////////////////////////////////////////////
  // JACOBIANS AND PREFACTORS 
  dsigma *= Jacob;

  // Changing units to zeptobarn = zb = 1e-45 cm^2
  dsigma *= (GeV2_to_cm2*1e45);
  
  myMC->xphys[0] = x1;
  myMC->xphys[1] = x2;
  myMC->xphys[2] = x3;
  myMC->xphys[3] = x4;
  myMC->xphys[4] = x5;
  myMC->xphys[5] = x6;
  myMC->xphys[6] = x7;
  myMC->xphys[7] = x8;
  myMC->xphys[8] = Enu;
  myMC->xphys[9] = dsigma;

  return 0;
}





long double hT_coh(int A, int Z, long double M, long double Enu, long double shat, long double Q2){


  long double hT = - 4 *  (Z*Z) * 4*M_PI*alpha_QED *(FF_WS(sqrt(Q2), A)*FF_WS(sqrt(Q2), A)) *
                            (shat*shat*M*M + 2*Q2*Enu*M*shat - 4*Q2*Enu*Enu*M*M );


  return hT;
}

long double hL_coh(int A, int Z, long double M, long double Enu, long double shat, long double Q2){


  long double hL = Q2 * (Z*Z) * 4 * M_PI * alpha_QED * (FF_WS(sqrt(Q2), A)*FF_WS(sqrt(Q2), A)) *
                            (4*Enu*M - shat)*(4*Enu*M - shat);


  return hL;
}
/*********************     						                   ***********************************************/
/*********************    Diff Cross section (nu GAMMA)  ***********************************************/
/*********************     				COHERENT               ***********************************************/

long double dsigma_dPS(int nu_alpha, int l1, int l2, int A, int Z, long double Enu, long double s, long double phi, long double theta, long double t, long double l, long double q, std::vector<long double> &BSM_params){


	// Define masses for l-(l1) and l+(l2)
	long double mj, mk, Vijk, Aijk;
	switch (l1)
		{
			case e_flag:
				mj = m_e;
				break;
			case mu_flag:
				mj = m_mu;
				break;
			case tau_flag:
				mj = m_tau;
				break;
		}

	switch (l2)
		{
			case e_flag:
				mk = m_e;
				break;
			case mu_flag:
				mk = m_mu;
				break;
			case tau_flag:
				mk = m_tau;
				break;
		}


	// Define the proper axial and vector coefficients for nu(nu_alpha), l-(l1) and l+(l2)
	if (nu_alpha == l1)
	{
		if (l1 == l2){Aijk = 0.5; Vijk = 0.5 + 2*sw2;}
		else if (l1 != l2) {Aijk = 1.0; Vijk = 1.0;}
		else {printf("Error! Flags for the leptons not well defined or not listed.");}
	}
	else if (nu_alpha != l1)
	{
		Aijk = -0.5; Vijk = -0.5 + 2*sw2;
	}
	else {printf("Error! Flags for the leptons not well defined or not listed.");}


  ///////////////////////////////////////////////////////////////////////////
  //               CHANGE OF VARIABLES TO CANCEL Z' PROPAGATOR
  ///////////////////////////////////////////////////////////////////////////
  

  // CHANGE OF VARIABLES!
  // l += t;

  if (BSM_params[0] == ZPRIME)
  {
    // W' and Z' extension under mu-tau gauge
    // long double gw=.5;// BSM W' coupling
    // long double gw=0;// BSM W' coupling
    // long double mwprime2=1e2*1e2; //mass squre of BSM W' boson in GeV*GeV
    long double mzprime = BSM_params[1]; //mass of BSM Z' boson in GeV
    long double gv=BSM_params[2];//  BSM Z' coupling to neutrino
    long double ga=BSM_params[3];// BSM Z' coupling  to charged lepton

    // if (nu_alpha == l1){
    //   Vijk=Vijk+4*gw*gw*mw2/mwprime2; 
    //   Aijk=Aijk+4*gw*gw*mw2/mwprime2; 
    // }
    if (l1 == l2){
      
      Vijk = Vijk - (SQR(gv)/sqrt(2.0)/Gf) * (1.0/(l-t-mzprime*mzprime)); // k^2 = l -t - q^2 
      Aijk = Aijk - (SQR(ga)/sqrt(2.0)/Gf) * (1.0/(l-t-mzprime*mzprime)); // k^2 = l - t - q^2

      // Vijk=Vijk-2*glz*gnuz*mW*mW/(l-t-q*q-mzprime*mzprime)/(1-sw2); // k^2 = l -t - q^2 
      // Vijk=Vijk-2*glz*gnuz*mW*mW/(l-t-mzprime*mzprime)/(1-sw2); // k^2 = l -t - q^2 
      // Aijk=Aijk-2*glz*gnuz*mW*mW/(l-t-q*q-mzprime*mzprime)/(1-sw2); // k^2 = l -t - q^2
      // Vijk=Vijk-2*glz*gnuz*mW*mW/(l/SQR(mzprime)-t/SQR(mzprime)-SQR(q/mzprime)-1.0)/(1-sw2)/SQR(mzprime); // k^2 = l -t - q^2 
      // Aijk=Aijk-2*glz*gnuz*mW*mW/(l/SQR(mzprime)-t/SQR(mzprime)-SQR(q/mzprime)-1.0)/(1-sw2)/SQR(mzprime); // k^2 = l -t - q^2
    }

    
  }
  
	// Effective mass *** NOTE SOME PLACES USE THE ALTERNATIVE Mjk = (Mj + Mk)/2 ***
	// long double ml1l2 = ml1 + ml2;

	long double sigma;

	long double mjkt = mj*mj + mk*mk;

  long double beta = sqrt(1.0 - 2.0*mjkt/l + pow(mj*mj - mk*mk,2)/l/l);

	sigma = //TA + TB
   128*pow(t,-4)*pow(-l + beta*l*cos(theta) - pow(mj,2) + pow(mk,2),-1)*pow(l + beta*l*cos(theta) - pow(mj,2) + pow(mk,2),-1)*
    (2*Aijk*l*t*Vijk*((l - t)*t*(-2*s + t)*(l - pow(mj,2) - pow(mk,2))*(pow(mj,2) - pow(mk,2)) + 
         beta*l*(l - t)*(-2*s + t)*cos(theta)*(2*pow(l,2) + t*(pow(mj,2) + pow(mk,2)) - 2*l*(t + pow(mj,2) + pow(mk,2))) + 
         beta*l*cos(phi)*(4*pow(l,2) - 2*l*(3*t + 2*pow(mj,2) + 2*pow(mk,2)) + t*(t + 4*pow(mj,2) + 4*pow(mk,2)))*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta)) + 
      pow(Aijk,2)*(pow(beta,2)*pow(l,2)*(-(s*(s - t)*(pow(mj,2) + pow(mk,2))*pow(t,2)) + pow(l,3)*pow(-2*s + t,2) - pow(l,2)*(2*t + pow(mj,2) + pow(mk,2))*pow(-2*s + t,2) + 
            l*t*(t*(-3*s*t + 3*pow(s,2) + pow(t,2)) + pow(mj,2)*pow(-2*s + t,2) + pow(mk,2)*pow(-2*s + t,2)))*pow(cos(theta),2) + 
         pow(t,2)*(pow(l,5) - pow(l,4)*(4*mj*mk + 2*t + pow(mj,2) + pow(mk,2)) + 
            pow(l,3)*(-(s*t) + 4*mk*pow(mj,3) - pow(mj,4) + 2*t*pow(mk,2) + 2*pow(mj,2)*(t + pow(mk,2)) + mj*(6*mk*t + 4*pow(mk,3)) - pow(mk,4) + pow(s,2) + pow(t,2)) + 
            pow(l,2)*(-4*mk*t*pow(mj,3) + pow(mj,6) + pow(mj,4)*(t - pow(mk,2)) + pow(mk,2)*(s*t + t*pow(mk,2) + pow(mk,4) - pow(s,2) - pow(t,2)) - 
               pow(mj,2)*(-(s*t) + 2*t*pow(mk,2) + pow(mk,4) + pow(s,2) + pow(t,2)) - 2*mj*mk*(-2*s*t + 2*t*pow(mk,2) + 2*pow(s,2) + pow(t,2))) - 
            s*(s - t)*(pow(mj,2) + pow(mk,2))*pow(pow(mj,2) - pow(mk,2),2) - l*(s*t + t*(pow(mj,2) + pow(mk,2)) - pow(s,2))*pow(pow(mj,2) - pow(mk,2),2)) - 
         4*s*(l - t)*(s - t)*pow(beta,2)*pow(l,3)*(l - t - pow(mj,2) - pow(mk,2))*pow(cos(phi),2)*pow(sin(theta),2) - 
         beta*l*cos(phi)*(t*(-2*s + t)*(pow(mj,2) - pow(mk,2))*(-2*l + t + 2*pow(mj,2) + 2*pow(mk,2)) + 
            2*beta*l*s*cos(theta)*(4*pow(l,2) + t*(t + 2*pow(mj,2) + 2*pow(mk,2)) - 2*l*(3*t + 2*pow(mj,2) + 2*pow(mk,2))))*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta) + 
         beta*l*t*cos(theta)*((pow(mj,2) - pow(mk,2))*(2*s*(s - t)*t*(pow(mj,2) + pow(mk,2)) + pow(l,2)*pow(-2*s + t,2) - 
               l*(4*s*(s - t)*pow(mj,2) + 4*s*(s - t)*pow(mk,2) + t*pow(-2*s + t,2))) + 
            beta*l*cos(phi)*(4*pow(l,2) + t*(t + 2*pow(mj,2) + 2*pow(mk,2)) - 2*l*(3*t + 2*pow(mj,2) + 2*pow(mk,2)))*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta))) + 
      pow(Vijk,2)*(pow(beta,2)*pow(l,2)*(-(s*(s - t)*(pow(mj,2) + pow(mk,2))*pow(t,2)) + pow(l,3)*pow(-2*s + t,2) - pow(l,2)*(2*t + pow(mj,2) + pow(mk,2))*pow(-2*s + t,2) + 
            l*t*(t*(-3*s*t + 3*pow(s,2) + pow(t,2)) + pow(mj,2)*pow(-2*s + t,2) + pow(mk,2)*pow(-2*s + t,2)))*pow(cos(theta),2) + 
         pow(t,2)*(pow(l,5) - pow(l,4)*(-4*mj*mk + 2*t + pow(mj,2) + pow(mk,2)) - 
            pow(l,3)*(s*t + 4*mk*pow(mj,3) + pow(mj,4) - 2*t*pow(mk,2) - 2*pow(mj,2)*(t + pow(mk,2)) + mj*(6*mk*t + 4*pow(mk,3)) + pow(mk,4) - pow(s,2) - pow(t,2)) + 
            pow(l,2)*(4*mk*t*pow(mj,3) + pow(mj,6) + pow(mj,4)*(t - pow(mk,2)) + pow(mk,2)*(s*t + t*pow(mk,2) + pow(mk,4) - pow(s,2) - pow(t,2)) - 
               pow(mj,2)*(-(s*t) + 2*t*pow(mk,2) + pow(mk,4) + pow(s,2) + pow(t,2)) + 2*mj*mk*(-2*s*t + 2*t*pow(mk,2) + 2*pow(s,2) + pow(t,2))) - 
            s*(s - t)*(pow(mj,2) + pow(mk,2))*pow(pow(mj,2) - pow(mk,2),2) - l*(s*t + t*(pow(mj,2) + pow(mk,2)) - pow(s,2))*pow(pow(mj,2) - pow(mk,2),2)) - 
         4*s*(l - t)*(s - t)*pow(beta,2)*pow(l,3)*(l - t - pow(mj,2) - pow(mk,2))*pow(cos(phi),2)*pow(sin(theta),2) - 
         beta*l*cos(phi)*(t*(-2*s + t)*(pow(mj,2) - pow(mk,2))*(-2*l + t + 2*pow(mj,2) + 2*pow(mk,2)) + 
            2*beta*l*s*cos(theta)*(4*pow(l,2) + t*(t + 2*pow(mj,2) + 2*pow(mk,2)) - 2*l*(3*t + 2*pow(mj,2) + 2*pow(mk,2))))*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta) + 
         beta*l*t*cos(theta)*((pow(mj,2) - pow(mk,2))*(2*s*(s - t)*t*(pow(mj,2) + pow(mk,2)) + pow(l,2)*pow(-2*s + t,2) - 
               l*(4*s*(s - t)*pow(mj,2) + 4*s*(s - t)*pow(mk,2) + t*pow(-2*s + t,2))) + 
            beta*l*cos(phi)*(4*pow(l,2) + t*(t + 2*pow(mj,2) + 2*pow(mk,2)) - 2*l*(3*t + 2*pow(mj,2) + 2*pow(mk,2)))*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta)))) + 
   32*pow(t,-4)*(pow(l - beta*l*cos(theta) + pow(mj,2) - pow(mk,2),-2)*((s - t)*t*(4*l*pow(mj,2) + t*(l - beta*l*cos(theta) + pow(mj,2) - pow(mk,2)))*pow(Aijk + Vijk,2)*
          ((l + s - t)*t*(l - pow(mj,2) + pow(mk,2)) - beta*l*cos(theta)*(2*l*s - l*t - s*t + pow(t,2)) + 2*beta*l*cos(phi)*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta)) + 
         2*pow(mj,2)*pow(Aijk + Vijk,2)*(-(beta*l*(2*l*s - l*t - s*t)*cos(theta)) + (l - s)*t*(l + pow(mj,2) - pow(mk,2)) + 2*beta*l*cos(phi)*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta))*
          ((l + s - t)*t*(l - pow(mj,2) + pow(mk,2)) - beta*l*cos(theta)*(2*l*s - l*t - s*t + pow(t,2)) + 2*beta*l*cos(phi)*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta)) - 
         (Aijk - Vijk)*(s*t*(Aijk - Vijk)*(4*l*pow(mj,2) + t*(l - beta*l*cos(theta) + pow(mj,2) - pow(mk,2)))*
             (beta*l*(2*l*s - l*t - s*t)*cos(theta) + (l - s)*t*(l - pow(mj,2) + pow(mk,2)) - 2*beta*l*cos(phi)*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta)) - 
            2*mj*(-2*l*mk*(l - t)*(Aijk + Vijk)*(4*l*pow(mj,2) + t*(-l + beta*l*cos(theta) - pow(mj,2) + pow(mk,2)))*pow(t,2) + 
               mj*(Aijk - Vijk)*(beta*l*(2*l*s - l*t - s*t)*cos(theta) + (l - s)*t*(l - pow(mj,2) + pow(mk,2)) - 2*beta*l*cos(phi)*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta))*
                ((l + s - t)*t*(l + pow(mj,2) - pow(mk,2)) + beta*l*cos(theta)*(2*l*s - l*t - s*t + pow(t,2)) - 2*beta*l*cos(phi)*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta))))) + 
      pow(l + beta*l*cos(theta) - pow(mj,2) + pow(mk,2),-2)*((s - t)*t*(4*l*pow(mk,2) + t*(l + beta*l*cos(theta) - pow(mj,2) + pow(mk,2)))*pow(Aijk - Vijk,2)*
          ((l + s - t)*t*(l + pow(mj,2) - pow(mk,2)) + beta*l*cos(theta)*(2*l*s - l*t - s*t + pow(t,2)) - 2*beta*l*cos(phi)*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta)) - 
         s*t*(4*l*pow(mk,2) + t*(l + beta*l*cos(theta) - pow(mj,2) + pow(mk,2)))*pow(Aijk + Vijk,2)*
          (-(beta*l*(2*l*s - l*t - s*t)*cos(theta)) + (l - s)*t*(l + pow(mj,2) - pow(mk,2)) + 2*beta*l*cos(phi)*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta)) + 
         2*mk*(mk*pow(Aijk + Vijk,2)*(-(beta*l*(2*l*s - l*t - s*t)*cos(theta)) + (l - s)*t*(l + pow(mj,2) - pow(mk,2)) + 2*beta*l*cos(phi)*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta))*
             ((l + s - t)*t*(l - pow(mj,2) + pow(mk,2)) - beta*l*cos(theta)*(2*l*s - l*t - s*t + pow(t,2)) + 2*beta*l*cos(phi)*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta)) + 
            (Aijk - Vijk)*(-2*l*mj*(l - t)*(Aijk + Vijk)*(4*l*pow(mk,2) - t*(l + beta*l*cos(theta) - pow(mj,2) + pow(mk,2)))*pow(t,2) + 
               mk*(Aijk - Vijk)*(beta*l*(2*l*s - l*t - s*t)*cos(theta) + (l - s)*t*(l - pow(mj,2) + pow(mk,2)) - 2*beta*l*cos(phi)*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta))*
                ((l + s - t)*t*(l + pow(mj,2) - pow(mk,2)) + beta*l*cos(theta)*(2*l*s - l*t - s*t + pow(t,2)) - 2*beta*l*cos(phi)*pow(l*s*(l - t)*(-s + t),0.5)*sin(theta))))));
  // nu gamma part!

	/* EXTRA FACTOR OF SIN(THETA) APPEARING IN theta INTEGRATION */
	sigma *= sin(theta) * alpha_QED*Gf*Gf/(16*pow(4.0*M_PI,3)*pow(s,2));

  /* EXTRA FACTOR OF BETA APPEARING IN l INTEGRATION */
	sigma *= beta;

  // HEAVISIDE FORM FACTOR 
  // sigma *= -2*log(s/(2.0*Enu*LAMBDA_QCD/pow(A,1.0/3.0)))/s;
  
  // SIMPLIFIED FORM FACTOR
  // sigma *= (-2*log(s/(2.0*Enu*0.143)) - log( 0.143*0.143 * (2)/( s*s/4/Enu/Enu + 0.143*0.143)   ))/s;

  // SIMPLIFIED FORM FACTOR
  // sigma *= SQR(1.0/(1.0+q*q/SQR(0.143))) * 2 /q/s;
  // sigma *= (2.0)* SQR( exp( - 0.5 * q*q / SQR( 0.197*sqrt(5)/1.2/pow(A,1.0/3.0)  )  ) )/q/s;

  // WOODS-SAXON FORM FACTOR
  sigma *= (2.0)*SQR(FF_WS(q, A))/q/s;


  // Changing units to zeptobarn = zb = 1e-45 cm^2
  sigma *= (GeV2_to_cm2*1e45);
  

	return sigma  * Z * Z * alpha_QED / M_PI;
}


/*********************                                   ***********************************************/
/*********************    Diff Cross section (nu GAMMA)  ***********************************************/
/*********************            COHERENT               ***********************************************/

long double dsigma_dPS_diff(int nu_alpha, int l1, int l2, int A, int Z, long double Enu, long double s, long double phi, long double theta, long double t, long double l, long double q, std::vector<long double> &BSM_params){


  // Define masses for l-(l1) and l+(l2)
  long double mj, mk, Vijk, Aijk;
  switch (l1)
    {
      case e_flag:
        mj = m_e;
        break;
      case mu_flag:
        mj = m_mu;
        break;
      case tau_flag:
        mj = m_tau;
        break;
    }

  switch (l2)
    {
      case e_flag:
        mk = m_e;
        break;
      case mu_flag:
        mk = m_mu;
        break;
      case tau_flag:
        mk = m_tau;
        break;
    }


  // Define the proper axial and vector coefficients for nu(nu_alpha), l-(l1) and l+(l2)
  if (nu_alpha == l1)
  {
    if (l1 == l2){Aijk = 0.5; Vijk = 0.5 + 2*sw2;}
    else if (l1 != l2) {Aijk = 1.0; Vijk = 1.0;}
    else {printf("Error! Flags for the leptons not well defined or not listed.");}
  }
  else if (nu_alpha != l1)
  {
    Aijk = -0.5; Vijk = -0.5 + 2*sw2;
  }
  else {printf("Error! Flags for the leptons not well defined or not listed.");}


  ///////////////////////////////////////////////////////////////////////////
  //               CHANGE OF VARIABLES TO CANCEL Z' PROPAGATOR
  ///////////////////////////////////////////////////////////////////////////

  if (BSM_params[0] == ZPRIME)
  {
    // W' and Z' extension under mu-tau gauge
    // long double gw=.5;// BSM W' coupling
    // long double gw=0;// BSM W' coupling
    // long double mwprime2=1e2*1e2; //mass squre of BSM W' boson in GeV*GeV
    long double mzprime = BSM_params[1]; //mass of BSM Z' boson in GeV
    long double gnuz = BSM_params[2];//  BSM Z' coupling to neutrino
    long double glz = BSM_params[3];// BSM Z' coupling  to charged lepton

    // if (nu_alpha == l1){
    //   Vijk=Vijk+4*gw*gw*mw2/mwprime2; 
    //   Aijk=Aijk+4*gw*gw*mw2/mwprime2; 
    // }
    if (l1 == l2){
      Vijk=Vijk - (SQR(glz)/sqrt(2.0)/Gf) * (1.0/(l-t-mzprime*mzprime)); // k^2 = l -t - q^2 
      Aijk = Aijk - (SQR(glz)/sqrt(2.0)/Gf) * (1.0/(l-t-mzprime*mzprime)); // k^2 = l - t - q^2

      // Vijk=Vijk-2*glz*gnuz*mW*mW/(l-t-q*q-mzprime*mzprime)/(1-sw2); // k^2 = l -t - q^2 
      // Aijk=Aijk-2*glz*gnuz*mW*mW/(l-t-q*q-mzprime*mzprime)/(1-sw2); // k^2 = l -t - q^2

      // Vijk=Vijk-2*glz*gnuz*mW*mW/(l/SQR(mzprime)-t/SQR(mzprime)-SQR(q/mzprime)-1.0)/(1-sw2)/SQR(mzprime); // k^2 = l -t - q^2 
      // Aijk=Aijk-2*glz*gnuz*mW*mW/(l/SQR(mzprime)-t/SQR(mzprime)-SQR(q/mzprime)-1.0)/(1-sw2)/SQR(mzprime); // k^2 = l -t - q^2
    }
  }

  // Effective mass *** NOTE SOME PLACES USE THE ALTERNATIVE Mjk = (Mj + Mk)/2 ***
  // long double ml1l2 = ml1 + ml2;

  long double sigma;

  long double mjkt = mj*mj + mk*mk;

  long double beta = sqrt(1.0 - 2.0*mjkt/l + pow(mj*mj - mk*mk,2)/l/l);

  sigma = //TA + TB
    (128*(2*Aijk*l*t*Vijk*((l - pow(mj,2) - pow(mk,2))*(pow(mj,2) - pow(mk,2))*(l - t)*t*(-2*s + t) + 
           beta*l*(l - t)*(-2*s + t)*(2*pow(l,2) + (pow(mj,2) + pow(mk,2))*t - 
              2*l*(pow(mj,2) + pow(mk,2) + t))*cos(theta) + 
           beta*l*sqrt(l*s*(l - t)*(-s + t))*
            (4*pow(l,2) + t*(4*pow(mj,2) + 4*pow(mk,2) + t) - 
              2*l*(2*pow(mj,2) + 2*pow(mk,2) + 3*t))*cos(phi)*sin(theta)) + 
        pow(Aijk,2)*(pow(t,2)*(pow(l,5) - 
              pow(pow(mj,2) - pow(mk,2),2)*(pow(mj,2) + pow(mk,2))*s*(s - t) - 
              pow(l,4)*(pow(mj,2) + 4*mj*mk + pow(mk,2) + 2*t) - 
              l*pow(pow(mj,2) - pow(mk,2),2)*(-pow(s,2) + (pow(mj,2) + pow(mk,2))*t + s*t) + 
              pow(l,3)*(-pow(mj,4) + 4*pow(mj,3)*mk - pow(mk,4) + pow(s,2) + 2*pow(mk,2)*t - 
                 s*t + pow(t,2) + 2*pow(mj,2)*(pow(mk,2) + t) + mj*(4*pow(mk,3) + 6*mk*t)) + 
              pow(l,2)*(pow(mj,6) - 4*pow(mj,3)*mk*t + pow(mj,4)*(-pow(mk,2) + t) + 
                 pow(mk,2)*(pow(mk,4) - pow(s,2) + pow(mk,2)*t + s*t - pow(t,2)) - 
                 2*mj*mk*(2*pow(s,2) + 2*pow(mk,2)*t - 2*s*t + pow(t,2)) - 
                 pow(mj,2)*(pow(mk,4) + pow(s,2) + 2*pow(mk,2)*t - s*t + pow(t,2)))) + 
           pow(beta,2)*pow(l,2)*(-((pow(mj,2) + pow(mk,2))*s*(s - t)*pow(t,2)) + 
              pow(l,3)*pow(-2*s + t,2) - 
              pow(l,2)*pow(-2*s + t,2)*(pow(mj,2) + pow(mk,2) + 2*t) + 
              l*t*(pow(mj,2)*pow(-2*s + t,2) + pow(mk,2)*pow(-2*s + t,2) + 
                 t*(3*pow(s,2) - 3*s*t + pow(t,2))))*pow(cos(theta),2) - 
           beta*l*sqrt(l*s*(l - t)*(-s + t))*cos(phi)*
            ((pow(mj,2) - pow(mk,2))*t*(-2*l + 2*pow(mj,2) + 2*pow(mk,2) + t)*(-2*s + t) + 
              2*beta*l*s*(4*pow(l,2) + t*(2*pow(mj,2) + 2*pow(mk,2) + t) - 
                 2*l*(2*pow(mj,2) + 2*pow(mk,2) + 3*t))*cos(theta))*sin(theta) - 
           4*pow(beta,2)*pow(l,3)*s*(l - t)*(l - pow(mj,2) - pow(mk,2) - t)*(s - t)*
            pow(cos(phi),2)*pow(sin(theta),2) + 
           beta*l*t*cos(theta)*((pow(mj,2) - pow(mk,2))*
               (2*(pow(mj,2) + pow(mk,2))*s*(s - t)*t + pow(l,2)*pow(-2*s + t,2) - 
                 l*(4*pow(mj,2)*s*(s - t) + 4*pow(mk,2)*s*(s - t) + t*pow(-2*s + t,2))) + 
              beta*l*sqrt(l*s*(l - t)*(-s + t))*
               (4*pow(l,2) + t*(2*pow(mj,2) + 2*pow(mk,2) + t) - 
                 2*l*(2*pow(mj,2) + 2*pow(mk,2) + 3*t))*cos(phi)*sin(theta))) + 
        pow(Vijk,2)*(pow(t,2)*(pow(l,5) - 
              pow(pow(mj,2) - pow(mk,2),2)*(pow(mj,2) + pow(mk,2))*s*(s - t) - 
              pow(l,4)*(pow(mj,2) - 4*mj*mk + pow(mk,2) + 2*t) - 
              l*pow(pow(mj,2) - pow(mk,2),2)*(-pow(s,2) + (pow(mj,2) + pow(mk,2))*t + s*t) - 
              pow(l,3)*(pow(mj,4) + 4*pow(mj,3)*mk + pow(mk,4) - pow(s,2) - 2*pow(mk,2)*t + 
                 s*t - pow(t,2) - 2*pow(mj,2)*(pow(mk,2) + t) + mj*(4*pow(mk,3) + 6*mk*t)) + 
              pow(l,2)*(pow(mj,6) + 4*pow(mj,3)*mk*t + pow(mj,4)*(-pow(mk,2) + t) + 
                 pow(mk,2)*(pow(mk,4) - pow(s,2) + pow(mk,2)*t + s*t - pow(t,2)) + 
                 2*mj*mk*(2*pow(s,2) + 2*pow(mk,2)*t - 2*s*t + pow(t,2)) - 
                 pow(mj,2)*(pow(mk,4) + pow(s,2) + 2*pow(mk,2)*t - s*t + pow(t,2)))) + 
           pow(beta,2)*pow(l,2)*(-((pow(mj,2) + pow(mk,2))*s*(s - t)*pow(t,2)) + 
              pow(l,3)*pow(-2*s + t,2) - 
              pow(l,2)*pow(-2*s + t,2)*(pow(mj,2) + pow(mk,2) + 2*t) + 
              l*t*(pow(mj,2)*pow(-2*s + t,2) + pow(mk,2)*pow(-2*s + t,2) + 
                 t*(3*pow(s,2) - 3*s*t + pow(t,2))))*pow(cos(theta),2) - 
           beta*l*sqrt(l*s*(l - t)*(-s + t))*cos(phi)*
            ((pow(mj,2) - pow(mk,2))*t*(-2*l + 2*pow(mj,2) + 2*pow(mk,2) + t)*(-2*s + t) + 
              2*beta*l*s*(4*pow(l,2) + t*(2*pow(mj,2) + 2*pow(mk,2) + t) - 
                 2*l*(2*pow(mj,2) + 2*pow(mk,2) + 3*t))*cos(theta))*sin(theta) - 
           4*pow(beta,2)*pow(l,3)*s*(l - t)*(l - pow(mj,2) - pow(mk,2) - t)*(s - t)*
            pow(cos(phi),2)*pow(sin(theta),2) + 
           beta*l*t*cos(theta)*((pow(mj,2) - pow(mk,2))*
               (2*(pow(mj,2) + pow(mk,2))*s*(s - t)*t + pow(l,2)*pow(-2*s + t,2) - 
                 l*(4*pow(mj,2)*s*(s - t) + 4*pow(mk,2)*s*(s - t) + t*pow(-2*s + t,2))) + 
              beta*l*sqrt(l*s*(l - t)*(-s + t))*
               (4*pow(l,2) + t*(2*pow(mj,2) + 2*pow(mk,2) + t) - 
                 2*l*(2*pow(mj,2) + 2*pow(mk,2) + 3*t))*cos(phi)*sin(theta)))))/
    (pow(t,4)*(-l - pow(mj,2) + pow(mk,2) + beta*l*cos(theta))*
      (l - pow(mj,2) + pow(mk,2) + beta*l*cos(theta))) + 
   (32*(((s - t)*t*pow(Aijk + Vijk,2)*(4*l*pow(mj,2) + 
              t*(l + pow(mj,2) - pow(mk,2) - beta*l*cos(theta)))*
            ((l - pow(mj,2) + pow(mk,2))*(l + s - t)*t - 
              beta*l*(2*l*s - l*t - s*t + pow(t,2))*cos(theta) + 
              2*beta*l*sqrt(l*s*(l - t)*(-s + t))*cos(phi)*sin(theta)) + 
           2*pow(mj,2)*pow(Aijk + Vijk,2)*
            ((l + pow(mj,2) - pow(mk,2))*(l - s)*t - beta*l*(2*l*s - l*t - s*t)*cos(theta) + 
              2*beta*l*sqrt(l*s*(l - t)*(-s + t))*cos(phi)*sin(theta))*
            ((l - pow(mj,2) + pow(mk,2))*(l + s - t)*t - 
              beta*l*(2*l*s - l*t - s*t + pow(t,2))*cos(theta) + 
              2*beta*l*sqrt(l*s*(l - t)*(-s + t))*cos(phi)*sin(theta)) - 
           (Aijk - Vijk)*(s*t*(Aijk - Vijk)*
               (4*l*pow(mj,2) + t*(l + pow(mj,2) - pow(mk,2) - beta*l*cos(theta)))*
               ((l - pow(mj,2) + pow(mk,2))*(l - s)*t + beta*l*(2*l*s - l*t - s*t)*cos(theta) - 
                 2*beta*l*sqrt(l*s*(l - t)*(-s + t))*cos(phi)*sin(theta)) - 
              2*mj*(-2*l*mk*(l - t)*pow(t,2)*(Aijk + Vijk)*
                  (4*l*pow(mj,2) + t*(-l - pow(mj,2) + pow(mk,2) + beta*l*cos(theta))) + 
                 mj*(Aijk - Vijk)*((l - pow(mj,2) + pow(mk,2))*(l - s)*t + 
                    beta*l*(2*l*s - l*t - s*t)*cos(theta) - 
                    2*beta*l*sqrt(l*s*(l - t)*(-s + t))*cos(phi)*sin(theta))*
                  ((l + pow(mj,2) - pow(mk,2))*(l + s - t)*t + 
                    beta*l*(2*l*s - l*t - s*t + pow(t,2))*cos(theta) - 
                    2*beta*l*sqrt(l*s*(l - t)*(-s + t))*cos(phi)*sin(theta)))))/
         pow(l + pow(mj,2) - pow(mk,2) - beta*l*cos(theta),2) + 
        ((s - t)*t*pow(Aijk - Vijk,2)*(4*l*pow(mk,2) + 
              t*(l - pow(mj,2) + pow(mk,2) + beta*l*cos(theta)))*
            ((l + pow(mj,2) - pow(mk,2))*(l + s - t)*t + 
              beta*l*(2*l*s - l*t - s*t + pow(t,2))*cos(theta) - 
              2*beta*l*sqrt(l*s*(l - t)*(-s + t))*cos(phi)*sin(theta)) - 
           s*t*pow(Aijk + Vijk,2)*(4*l*pow(mk,2) + 
              t*(l - pow(mj,2) + pow(mk,2) + beta*l*cos(theta)))*
            ((l + pow(mj,2) - pow(mk,2))*(l - s)*t - beta*l*(2*l*s - l*t - s*t)*cos(theta) + 
              2*beta*l*sqrt(l*s*(l - t)*(-s + t))*cos(phi)*sin(theta)) + 
           2*mk*(mk*pow(Aijk + Vijk,2)*((l + pow(mj,2) - pow(mk,2))*(l - s)*t - 
                 beta*l*(2*l*s - l*t - s*t)*cos(theta) + 
                 2*beta*l*sqrt(l*s*(l - t)*(-s + t))*cos(phi)*sin(theta))*
               ((l - pow(mj,2) + pow(mk,2))*(l + s - t)*t - 
                 beta*l*(2*l*s - l*t - s*t + pow(t,2))*cos(theta) + 
                 2*beta*l*sqrt(l*s*(l - t)*(-s + t))*cos(phi)*sin(theta)) + 
              (Aijk - Vijk)*(-2*l*mj*(l - t)*pow(t,2)*(Aijk + Vijk)*
                  (4*l*pow(mk,2) - t*(l - pow(mj,2) + pow(mk,2) + beta*l*cos(theta))) + 
                 mk*(Aijk - Vijk)*((l - pow(mj,2) + pow(mk,2))*(l - s)*t + 
                    beta*l*(2*l*s - l*t - s*t)*cos(theta) - 
                    2*beta*l*sqrt(l*s*(l - t)*(-s + t))*cos(phi)*sin(theta))*
                  ((l + pow(mj,2) - pow(mk,2))*(l + s - t)*t + 
                    beta*l*(2*l*s - l*t - s*t + pow(t,2))*cos(theta) - 
                    2*beta*l*sqrt(l*s*(l - t)*(-s + t))*cos(phi)*sin(theta)))))/
         pow(l - pow(mj,2) + pow(mk,2) + beta*l*cos(theta),2)))/pow(t,4);
  // nu gamma part!

  /* EXTRA FACTOR OF SIN(THETA) APPEARING IN theta INTEGRATION */
  sigma *= sin(theta) * alpha_QED*Gf*Gf/(16*pow(4.0*M_PI,3)*pow(s,2));

  /* EXTRA FACTOR OF BETA APPEARING IN l INTEGRATION */
  sigma *= beta;

  // DIRAC DIPOLE FORM FACTOR
  sigma *= SQR(F_diffractive(q*q))*2/q/s;

  // Changing units to zeptobarn = zb = 1e-45 cm^2
  sigma *= (GeV2_to_cm2*1e45);

  return sigma  * alpha_QED / M_PI ;
}

/***********************     						          ***********************************************/
/***********************        FORM FACTORS      ***********************************************/
/***********************     						          ***********************************************/


/***********************        COHERENT      ***********************************************/
long double FF_WS(long double q, long double A){

  long double r0 = 1.07/0.197 * pow(A,1.0/3.0);
  long double a = 0.57/0.197;

  return 3.0*M_PI*q*a*( M_PI*q*a*(1.0/tanh(M_PI*q*a))*sin(q*r0) - q*r0*cos(q*r0) )/ (q*q*q*sinh(M_PI*q*a)*(r0*r0+M_PI*M_PI*a*a)*r0);
}

long double u1_WS(long double x1, long double A){

  long double beta2 = SQR(sqrt(5)*0.197/1.2/pow(A,1.0/3.0));
  return exp(- x1/beta2) * beta2;
}

long double x1_WS(long double u1, long double A){

  long double beta2 = SQR(sqrt(5)*0.197/1.2/pow(A,1.0/3.0));
  return - beta2 * log( u1/beta2 );
}

/***********************        DIFFRACTIVE      ***********************************************/

long double F_1(long double q2){

  long double G_dip = pow( (1.0+ q2/0.71), -2.0); // 0.71 GeV^2
  long double tau_param = q2/4.0/SQR(M_avg);

  return (G_dip + tau_param*qsi*G_dip) / (1.0 + tau_param);
}


long double H1_n(long double q){
  long double tau_param = -SQR(q)/4.0/SQR(M_avg); // Recall q = sqrt(-q^2) = |Q^2|, so tau = -|Q|^2/4M^2
  return -tau_param/(1.0 - tau_param)*SQR(MAG_MOM_N*G_dip(q));
}
long double H2_n(long double q){
  long double tau_param = -SQR(q)/4.0/SQR(M_avg); // Recall q = sqrt(-q^2) = |Q^2|, so tau = -|Q|^2/4M^2
  return SQR(MAG_MOM_N*G_dip(q));
}
long double H1_p(long double q){
  long double tau_param = -SQR(q)/4.0/SQR(M_avg); // Recall q = sqrt(-q^2) = |Q^2|, so tau = -|Q|^2/4M^2
  return (SQR(G_dip(q))- tau_param*SQR(MAG_MOM_P*G_dip(q)))/(1.0 - tau_param); 
}
long double H2_p(long double q){
  long double tau_param = -SQR(q)/4.0/SQR(M_avg); // Recall q = sqrt(-q^2) = |Q^2|, so tau = -|Q|^2/4M^2
  return SQR(MAG_MOM_P*G_dip(q));
}


/***********************        OLD ONES      ***********************************************/

long double F(long double q2, long double A){

  // HEAVISIDE function that depends on the NUCLEAR MASS via A Ar_18^40)
  long double q2MAX = SQR(LAMBDA_QCD/pow(A,1.0/3.0));
  return heaviside(q2MAX - q2);
}

long double F_diffractive(long double q2){
  long double G_d = pow( (1.0+ q2/0.71), -2.0); // 0.71 GeV^2
  long double tau_param = q2/4.0/SQR(M_avg);
  return (G_d + tau_param*qsi*G_d) / (1.0 + tau_param);
}

/***********************     PAULI BLOCKING   ***********************************************/
long double F_PAULI_BLOCK(long double x1){
  
  long double qvec = sqrt( x1*(1+x1/4.0/SQR(M_avg))  );

  return (1.5*(qvec/2.0/P_FERMI) - 0.5*(qvec/2.0/P_FERMI)*(qvec/2.0/P_FERMI)*(qvec/2.0/P_FERMI) )*(qvec <= 2*P_FERMI) 
                + 1*(qvec > 2*P_FERMI);
      
}

