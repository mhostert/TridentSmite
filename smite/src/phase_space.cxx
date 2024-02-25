#include "constants.h"
#include "cross_sections.h"
#include "integrator.h"
#include "phase_space.h"

std::vector<long double> phase_space_flux(const cubareal xx[], void *MC){

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

	  myMC->V2     = myMC->Vijk*myMC->Vijk;
	  myMC->A2     = myMC->Aijk*myMC->Aijk;
	  myMC->VA     = myMC->Vijk*myMC->Aijk;
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

	  myMC->V2     = myMC->Vijk*myMC->Vijk * prop;
	  myMC->A2     = myMC->Aijk*myMC->Aijk * prop;
	  myMC->VA     = myMC->Vijk*myMC->Aijk * prop;
	  break;

	case (BSMonly):

	  u6_l = 1.0/(2.0*m6_u + mzprime*mzprime);
	  u6_u = 1.0/(2.0*m6_l + mzprime*mzprime);
	  u6   = u6_s*(u6_u - u6_l) + u6_l;
	  m6   = -mzprime*mzprime/2.0 + 1.0/2.0/u6;

	  prop = -u6;
	  jacob_u6 = 1.0/(u6*u6)/2.0;

	  myMC->V2     = myMC->Vijk*myMC->Vijk * prop*prop;
	  myMC->A2     = myMC->Aijk*myMC->Aijk * prop*prop;
	  myMC->VA     = myMC->Vijk*myMC->Aijk * prop*prop;
	  break;

	case (SMandBSM):
	 
	  u6_l = m6_l;
	  u6_u = m6_u;
	  u6   = u6_s*(u6_u - u6_l) + u6_l;
	  m6 = u6;
	  jacob_u6 = 1.0;

	  prop = 1.0/(2.0*m6 + mzprime*mzprime);

	  myMC->V2     = (myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop) ;
	  myMC->A2     = (myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop) ;
	  myMC->VA     = (myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop);
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
	// long double C4 = (E4*q0 - x4)/(p4vec*qvec);
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


	std::vector<long double> result(10);

	result[0] = x1;
	result[1] = x2;
	result[2] = x3;
	result[3] = x4;
	result[4] = x5;
	result[5] = x6;
	result[6] = x7;
	result[7] = x8;
	result[8] = Enu;
	result[9] = Jacob;

	return result;	
}



std::vector<long double> phase_space_Efixed(const cubareal xx[], void *MC){

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

	  myMC->V2     = myMC->Vijk*myMC->Vijk;
	  myMC->A2     = myMC->Aijk*myMC->Aijk;
	  myMC->VA     = myMC->Vijk*myMC->Aijk;
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

	  myMC->V2     = myMC->Vijk*myMC->Vijk * prop;
	  myMC->A2     = myMC->Aijk*myMC->Aijk * prop;
	  myMC->VA     = myMC->Vijk*myMC->Aijk * prop;
	  break;

	case (BSMonly):

	  u6_l = 1.0/(2.0*m6_u + mzprime*mzprime);
	  u6_u = 1.0/(2.0*m6_l + mzprime*mzprime);
	  u6   = u6_s*(u6_u - u6_l) + u6_l;
	  m6   = -mzprime*mzprime/2.0 + 1.0/2.0/u6;

	  prop = -u6;
	  jacob_u6 = 1.0/(u6*u6)/2.0;

	  myMC->V2     = myMC->Vijk*myMC->Vijk * prop*prop;
	  myMC->A2     = myMC->Aijk*myMC->Aijk * prop*prop;
	  myMC->VA     = myMC->Vijk*myMC->Aijk * prop*prop;
	  break;

	case (SMandBSM):
	 
	  u6_l = m6_l;
	  u6_u = m6_u;
	  u6   = u6_s*(u6_u - u6_l) + u6_l;
	  m6 = u6;
	  jacob_u6 = 1.0;

	  prop = 1.0/(2.0*m6 + mzprime*mzprime);

	  myMC->V2     = (myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop) ;
	  myMC->A2     = (myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop) ;
	  myMC->VA     = (myMC->Aijk + myMC->gprimeA*myMC->gprimeA/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop)*(myMC->Vijk + myMC->gprimeV*myMC->gprimeV/2.0/sqrt(2.0)/Gf * myMC->CHARGE * prop);
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
	// long double C4 = (E4*q0 - x4)/(p4vec*qvec);
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


	std::vector<long double> result(10);

	result[0] = x1;
	result[1] = x2;
	result[2] = x3;
	result[3] = x4;
	result[4] = x5;
	result[5] = x6;
	result[6] = x7;
	result[7] = x8;
	result[8] = Jacob;

	return result;	
}





def S_to_LAB_params():

  long double Q2   = x[0];      // -q^2 x1
  long double dot_q1   = x[1];      // q.k1 dot_q1
  long double Enu   = x[8];     // Enu in lab

  /************************************************************************************************/
  long double mN = Mnarg; //Nuclear mass.
    
  long double q0 = -Q2/(2*mN);
  long double qP = sqrt( Q2 + Q2*Q2/(4*mN*mN) );
  long double cosq = (Enu*q0 - dot_q1)/(qP*Enu);
  long double sinq = sqrt(1.-cosq*cosq);

  long double alpha = atan2(qP*sinq,(Enu+qP*cosq));
  long double sina = sin(alpha);
  long double cosa = cos(alpha);
  long double beta = (sina*qP*sinq + cosa*(Enu+qP*cosq))/(Enu+q0);
  long double gamma = 1.0/sqrt(1.-beta*beta);


  long double delta = atan2( -sina,(gamma*(cosa-beta))  );
  long double sind = sin(delta);
  long double cosd = cos(delta);

  std::vector<long double> params(6);


  params[0] = gamma;
  params[1] = beta;
  params[2] = sina;
  params[3] = cosa;
  params[4] = sind;
  params[5] = cosd;

  return params;


def P_Sframe(self, samples):

	Q2   = x[0];      // -q^2 x1
	dot_q1   = x[1];      // q.k1 x2
	dot_qplus    = x[2]; // q.p_plus x3
	dot_qminus   = x[3]; // q.p_minus x4
	dot_1plus    = x[4]; // k1.p_plus x5
	dot_1minus   = x[5]; // k1.p_minus x6
	phi_minus = x[6]; // x7 random number (0,2pi) for the special-frame azimuthal angle (taken to be the l_minus azimuthal angle in that frame).
	// long double phi = x[7]; // x8 random number (0,2pi) for the trivial overall azimuthal angle
	// long double Enu   = x[8];     // Enu in lab

	/************************************************************************************************/
	m_minus = mjarg;
	m_plus  = mkarg;

	shat = 2*dot_q1 - Q2; // (k1 + q) COM energy
	dot_12 = dot_q1 - dot_1plus - dot_1minus; // k1.p2
	dot_q2 = dot_q1 - dot_qplus - dot_qminus - Q2; // q.p2
	Eplus = (dot_1plus + dot_qplus)/sqrt(shat);
	Eminus = (dot_1minus + dot_qminus)/sqrt(shat);
	Ep2 = (dot_12+dot_q2)/sqrt(shat);

	Pplus = sqrt(Eplus*Eplus - m_plus*m_plus);
	Pminus = sqrt(Eminus*Eminus - m_minus*m_minus);
	Pp2 = Ep2;

	Eq = (dot_q1-Q2)/(sqrt(shat)); // energy of q in S frame
	p1 = (dot_q1)/(sqrt(shat)); // momentum of k1 and -q in S frame

	cosplus = (p1*dot_qplus - Eq*dot_1plus)/(p1*Pplus*sqrt(shat));
	sinplus = sqrt(1-cosplus*cosplus);
	cosminus = (p1*dot_qminus - Eq*dot_1minus)/(p1*Pminus*sqrt(shat));
	sinminus = sqrt(1-cosminus*cosminus);

	cosp2 = (p1*dot_q2 - Eq*dot_12)/(p1*Pp2*sqrt(shat));
	sinp2 = sqrt(1-cosp2*cosp2);

	// We have to work a bit harder to get the azmithual angle of phi_plus from phi_minus.
	// The z_ij = dot_ij 
	dot_pm = (-1)*(dot_q1 - dot_qplus - dot_qminus - dot_1plus - dot_1minus + 0.5*(m_minus*m_minus + m_plus*m_plus - Q2));
	dot_m2 = (dot_q1 - dot_qplus - dot_1plus - 0.5*(m_minus*m_minus - m_plus*m_plus + Q2));
	/************************************************************************************************/

	// The sign of p2 does not matter it seems.
	SIGN_PM = +1;
	SIGN_M2 = -1;
	phi_plus = phi_minus + SIGN_PM * acos((Eplus*Eminus - dot_pm - Pplus*Pminus*cosplus*cosminus)/(Pplus*Pminus*sinplus*sinminus));
	phi_p2 = phi_minus + SIGN_M2 * acos((Ep2*Eminus - dot_m2 - Pp2*Pminus*cosp2*cosminus)/(Pp2*Pminus*sinp2*sinminus));

	std::vector<long double> P2(4), P3(4), P4(4);
	std::vector<std::vector<long double>> vector_of_Ps;

	// E
	P2[0] = Ep2;
	// Px
	P2[1] = Pp2*sinp2*cos(phi_p2);
	// Py
	P2[2] = Pp2*sinp2*sin(phi_p2);
	// Pz
	P2[3] = Pp2*cosp2;
		
	// E
	P3[0] = Eplus;
	// Px
	P3[1] = Pplus*sinplus*cos(phi_plus);
	// Py  
	P3[2] = Pplus*sinplus*sin(phi_plus);
	// Pz
	P3[3] = Pplus*cosplus;

	// E
	P4[0] = Eminus;
	// Px
	P4[1] = Pminus*sinminus*cos(phi_minus);
	// Py
	P4[2] = Pminus*sinminus*sin(phi_minus);
	// Pz
	P4[3] = Pminus*cosminus;

	vector_of_Ps.push_back(P2);
	vector_of_Ps.push_back(P3);
	vector_of_Ps.push_back(P4);

	return vector_of_Ps;
