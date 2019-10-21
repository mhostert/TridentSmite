#include "observables.h"
  

/*********************                                                ***********************************************/
/*********************    4 MOMENTUM OF l3 or l^+                     ***********************************************/
/*********************                                                ***********************************************/
std::vector<long double> S_to_LAB(const std::vector<long double> &x, const std::vector<long double> &params){

  std::vector<long double> x_L;

  long double gamma = params[0];
  long double beta = params[1];
  long double sina = params[2];
  long double cosa = params[3];
  long double sind = params[4];
  long double cosd = params[5];

  x_L.push_back( gamma*x[0] - gamma*beta*(sind*x[1] - cosd*x[3])    );
  x_L.push_back( gamma*beta*sina*x[0] + (cosa*cosd -gamma*sina*sind)*x[1] + (cosa*sind+sina*cosd*gamma)*x[3]  );
  x_L.push_back( x[2] );
  x_L.push_back( gamma*beta*cosa*x[0] - (sina*cosd + gamma*cosa*sind)*x[1] - (sina*sind-gamma*cosa*cosd)*x[3]  );

  return x_L;
}

std::vector<long double> S_to_LAB_params(const long double mjarg, const long double mkarg, const long double Mnarg, const std::vector<long double> &x){

  long double x1   = x[0];      // -q^2 x1
  long double x2   = x[1];      // q.k1 x2
  long double b_plus    = x[2]; // q.p_plus x3
  long double b_minus   = x[3]; // q.p_minus x4
  long double a_plus    = x[4]; // k1.p_plus x5
  long double a_minus   = x[5]; // k1.p_minus x6
  long double phi_minus = x[6]; // a random number (0,2pi) for the special-frame azimuthal angle (taken to be the l_minus azimuthal angle in that frame).
  long double phi = x[7]; // a random number (0,2pi) for the trivial overall azimuthal angle
  long double Enu   = x[8];     // Enu in lab

  /************************************************************************************************/
  long double mN = Mnarg; //Nuclear mass.
    
  long double m_minus = mjarg;
  long double m_plus  = mkarg;

  long double q0 = -x1/(2*mN);
  long double qP = sqrt( x1 + x1*x1/(4*mN*mN) );
  long double cosq = (Enu*q0 - x2)/(qP*Enu);
  long double sinq = sqrt(1.-cosq*cosq);
  if (abs(sinq)>1.0 || abs(cosq)>1.0)
  {
    std::cout<<"ERROR! Sinq or cosq."<<std::endl;
  }

  long double alpha = atan2(qP*sinq,(Enu+qP*cosq));
  long double sina = sin(alpha);
  long double cosa = cos(alpha);
  long double beta = (sina*qP*sinq + cosa*(Enu+qP*cosq))/(Enu+q0);
  long double gamma = 1.0/sqrt(1.-beta*beta);
  if (abs(sina)>1.0 || abs(cosa)>1.0)
  {
    std::cout<<"ERROR! SIN or COS of a."<<std::endl;
  }

  long double delta = atan2( -sina,(gamma*(cosa-beta))  );
  long double sind = sin(delta);
  long double cosd = cos(delta);
  if (abs(sind)>1.0 || abs(cosd)>1.0)
  {
    std::cout<<"ERROR! SIN or COS of d."<<std::endl;
  }

  long double shat = 2*x2 - x1;

  long double a2 = x2 - a_plus - a_minus;
  long double b2 = x2 - b_plus - b_minus - x1;

  long double Eplus = (a_plus + b_plus)/sqrt(shat);
  long double Eminus = (a_minus + b_minus)/sqrt(shat);
  long double Ep2 = (a2+b2)/sqrt(shat);

  long double Pplus = sqrt(Eplus*Eplus - m_plus*m_plus);
  long double Pminus = sqrt(Eminus*Eminus - m_minus*m_minus);
  long double Pp2 = sqrt(Ep2*Ep2);

  long double eta0 = (shat - x1)/(2*sqrt(shat));
  long double epsilon1 = (shat + x1)/(2*sqrt(shat));
  long double cosplus = (epsilon1*b_plus - eta0*a_plus)/(epsilon1*Pplus*sqrt(shat));
  long double sinplus = sqrt(1.0-cosplus*cosplus);
  long double cosminus = (epsilon1*b_minus - eta0*a_minus)/(epsilon1*Pminus*sqrt(shat));
  long double sinminus = sqrt(1-cosminus*cosminus);
  
  // if (abs(cosplus)>1.0 || abs(sinplus)>1.0)
  // {
  //   std::cout<<"ERROR! SIN or COS of plus."<<std::endl;
  //   std::cout<<cosplus<<std::endl;
  // }

  // if (abs(cosminus)>1.0 || abs(sinminus)>1.0)
  // {
  //   std::cout<<"ERROR! SIN or COS of minus."<<std::endl;
  // }

  long double cosp2 = (epsilon1*b2 - eta0*a2)/(epsilon1*Pp2*sqrt(shat));
  long double sinp2 = sqrt(1-cosp2*cosp2);
  
  if (abs(sinp2)>1.0 || abs(cosp2)>1.0)
  {
    std::cout<<"ERROR! SIN or COS of p2."<<std::endl;
  }

  // We have to work a bit harder to get the azmithual angle of phi_plus from phi_minus.
  // The z_inv here is the Pplus.Pminus invariant.
  long double z_pm = (-1)*(x2 - b_plus - b_minus - a_plus - a_minus + 0.5*(m_minus*m_minus + m_plus*m_plus - x1));
  long double z_m2 = (x2 - b_plus - a_plus - 0.5*(m_minus*m_minus - m_plus*m_plus + x1));
  // double z_inv = x2 - b_plus - b_minus - a_plus - a_minus +0.5*(m_minus*m_minus + m_plus*m_plus - x1);
  // long double z_inv = (-1)*(m6 - b_plus - b_minus +0.5*(m_minus*m_minus + m_plus*m_plus - x1));
  /************************************************************************************************/
  std::vector<long double> params;


  params.push_back(gamma);
  params.push_back(beta);
  params.push_back(sina);
  params.push_back(cosa);
  params.push_back(sind);
  params.push_back(cosd);

  return params;
}


std::vector<long double> P_plus_Sframe(const long double mjarg, const long double mkarg, const long double Mnarg, const std::vector<long double> &x){

  long double x1   = x[0];      // -q^2 x1
  long double x2   = x[1];      // q.k1 x2
  long double b_plus    = x[2]; // q.p_plus x3
  long double b_minus   = x[3]; // q.p_minus x4
  long double a_plus    = x[4]; // k1.p_plus x5
  long double a_minus   = x[5]; // k1.p_minus x6
  long double phi_minus = x[6]; // a random number (0,2pi) for the special-frame azimuthal angle (taken to be the l_minus azimuthal angle in that frame).
  long double phi = x[7]; // a random number (0,2pi) for the trivial overall azimuthal angle
  long double Enu   = x[8];     // Enu in lab

  /************************************************************************************************/
  long double mN = Mnarg; //Nuclear mass.
  long double m_minus = mjarg;
  long double m_plus  = mkarg;

  long double shat = 2*x2 - x1;
  long double a2 = x2 - a_plus - a_minus;
  long double b2 = x2 - b_plus - b_minus - x1;
  long double Eplus = (a_plus + b_plus)/sqrt(shat);
  long double Eminus = (a_minus + b_minus)/sqrt(shat);
  long double Ep2 = (a2+b2)/sqrt(shat);

  long double Pplus = sqrt(Eplus*Eplus - m_plus*m_plus);
  long double Pminus = sqrt(Eminus*Eminus - m_minus*m_minus);
  long double Pp2 = sqrt(Ep2*Ep2);

  long double eta0 = (x2-x1)/(sqrt(shat));
  long double epsilon1 = (x2)/(sqrt(shat));
  long double cosplus = (epsilon1*b_plus - eta0*a_plus)/(epsilon1*Pplus*sqrt(shat));
  long double sinplus = sqrt(1-cosplus*cosplus);
  long double cosminus = (epsilon1*b_minus - eta0*a_minus)/(epsilon1*Pminus*sqrt(shat));
  long double sinminus = sqrt(1-cosminus*cosminus);
  
  // if (abs(cosplus)>1.0 || abs(sinplus)>1.0)
  // {
  //   std::cout<<"ERROR! SIN or COS of plus."<<std::endl;
  // }

  // if (abs(cosminus)>1.0 || abs(sinminus)>1.0)
  // {
  //   std::cout<<"ERROR! SIN or COS of minus."<<std::endl;
  // }

  long double cosp2 = (epsilon1*b2 - eta0*a2)/(epsilon1*Pp2*sqrt(shat));
  long double sinp2 = sqrt(1-cosp2*cosp2);
  
  if (abs(sinp2)>1.0 || abs(cosp2)>1.0)
  {
    std::cout<<"ERROR! SIN or COS of p2."<<std::endl;
  }

  // We have to work a bit harder to get the azmithual angle of phi_plus from phi_minus.
  // The z_inv here is the Pplus.Pminus invariant.
  long double z_pm = (-1)*(x2 - b_plus - b_minus - a_plus - a_minus + 0.5*(m_minus*m_minus + m_plus*m_plus - x1));
  long double z_m2 = (x2 - b_plus - a_plus - 0.5*(m_minus*m_minus - m_plus*m_plus + x1));
  // double z_inv = x2 - b_plus - b_minus - a_plus - a_minus +0.5*(m_minus*m_minus + m_plus*m_plus - x1);
  // long double z_inv = (-1)*(m6 - b_plus - b_minus +0.5*(m_minus*m_minus + m_plus*m_plus - x1));
  /************************************************************************************************/

  long double SIGN_PM = +1;
  long double SIGN_M2 = -1;
  long double phi_plus = phi_minus + SIGN_PM * acos((Eplus*Eminus - z_pm - Pplus*Pminus*cosplus*cosminus)/(Pplus*Pminus*sinplus*sinminus));
  long double phi_p2 = phi_minus + SIGN_M2 * acos((Ep2*Eminus - z_m2 - Pp2*Pminus*cosp2*cosminus)/(Pp2*Pminus*sinp2*sinminus));
  
  std::vector<long double> P2, P3;

  // E
  P2.push_back( Ep2  );
  // Px
  P2.push_back( Pp2*sinp2*cos(phi_p2)  );
  // Py
  P2.push_back( Pp2*sinp2*sin(phi_p2) );
  // Pz
  P2.push_back( Pp2*cosp2 );
    
  // E
  P3.push_back( Eplus  );
  // Px
  P3.push_back( Pplus*sinplus*cos(phi_plus)  );
  // Py  
  P3.push_back( Pplus*sinplus*sin(phi_plus) );
  // Pz
  P3.push_back( Pplus*cosplus );



  return P3;
}


std::vector<long double> P_minus_Sframe(const long double mjarg, const long double mkarg, const long double Mnarg, const std::vector<long double> &x){

  long double x1   = x[0];      // -q^2 x1
  long double x2   = x[1];      // q.k1 x2
  long double b_plus    = x[2]; // q.p_plus x3
  long double b_minus   = x[3]; // q.p_minus x4
  long double a_plus    = x[4]; // k1.p_plus x5
  long double a_minus   = x[5]; // k1.p_minus x6
  long double phi_minus = x[6]; // a random number (0,2pi) for the special-frame azimuthal angle (taken to be the l_minus azimuthal angle in that frame).
  long double phi = x[7]; // a random number (0,2pi) for the trivial overall azimuthal angle
  long double Enu   = x[8];     // Enu in lab

  /************************************************************************************************/
  long double mN = Mnarg; //Nuclear mass.
  long double m_minus = mjarg;
  long double m_plus  = mkarg;

  long double shat = 2*x2 - x1;
  long double a2 = x2 - a_plus - a_minus;
  long double b2 = x2 - b_plus - b_minus - x1;
  long double Eplus = (a_plus + b_plus)/sqrt(shat);
  long double Eminus = (a_minus + b_minus)/sqrt(shat);
  long double Ep2 = (a2+b2)/sqrt(shat);


  long double Pplus = sqrt(Eplus*Eplus - m_plus*m_plus);
  long double Pminus = sqrt(Eminus*Eminus - m_minus*m_minus);
  long double Pp2 = sqrt(Ep2*Ep2);

  long double eta0 = (x2-x1)/(sqrt(shat));
  long double epsilon1 = (x2)/(sqrt(shat));
  long double cosplus = (epsilon1*b_plus - eta0*a_plus)/(epsilon1*Pplus*sqrt(shat));
  long double sinplus = sqrt(1.-cosplus*cosplus);
  long double cosminus = (epsilon1*b_minus - eta0*a_minus)/(epsilon1*Pminus*sqrt(shat));
  long double sinminus = sqrt(1.-cosminus*cosminus);
  
  // if (abs(cosplus)>1.0 || abs(sinplus)>1.0)
  // {
  //   std::cout<<"ERROR! SIN or COS of plus."<<std::endl;
  // }

  // if (abs(cosminus)>1.0 || abs(sinminus)>1.0)
  // {
  //   std::cout<<"ERROR! SIN or COS of minus."<<std::endl;
  // }

  long double cosp2 = (epsilon1*b2 - eta0*a2)/(epsilon1*Pp2*sqrt(shat));
  long double sinp2 = sqrt(1-cosp2*cosp2);
  
  if (abs(sinp2)>1.0 || abs(cosp2)>1.0)
  {
    std::cout<<"ERROR! SIN or COS of p2."<<std::endl;
  }

  std::vector<long double> P4;

  // E
  P4.push_back( Eminus  );
  // Px
  P4.push_back( Pminus*sinminus*cos(phi_minus)  );
  // Py
  P4.push_back( Pminus*sinminus*sin(phi_minus) );
  // Pz
  P4.push_back( Pminus*cosminus );

  return P4;
}


////////////////////////////////////////////////////////////////////////////////////////
//     OTHER 4 MOMENTA!!! The outgoing neutrino P2 is in the P_plus_Sframe function
////////////////////////////////////////////////////////////////////////////////////////
  // std::vector<long double> P4,P4_L,P1_L,P1,q,PN_L,q_L;

////////////////////////////////////////////////////////////////////////////////////////
//     OTHER 4 MOMENTA!!! The outgoing neutrino P2 is in the P_plus_Sframe function
////////////////////////////////////////////////////////////////////////////////////////
  // std::vector<long double> P4,P4_L,P1_L,P1,q,PN_L,q_L;

  // // E 
  // P4.push_back( Eminus  );
  // // Px
  // P4.push_back( Pminus*sinminus*cos(phi_minus)  );
  // // Py
  // P4.push_back( Pminus*sinminus*sin(phi_minus) );
  // // Pz
  // P4.push_back( Pminus*cosminus );
     

  // P4_L = S_to_LAB(P4, params);

  // // NoP4 P4e need to boost it all back into the lab frame.
  // P1_L.push_back(Enu);
  // P1_L.push_back(0.0);
  // P1_L.push_back(0.0);
  // P1_L.push_back(Enu);

  // P1.push_back(epsilon1);
  // P1.push_back(0.0);
  // P1.push_back(0.0);
  // P1.push_back(epsilon1);



  // // NoP4 P4e need to boost it all back into the lab frame.
  // PN_L.push_back(mN);
  // PN_L.push_back(0.0);
  // PN_L.push_back(0.0);
  // PN_L.push_back(0.0);

  // q_L.push_back(q0);
  // q_L.push_back(qP*sinq);
  // q_L.push_back(0);
  // q_L.push_back(qP*cosq);

  // q.push_back(eta0);
  // q.push_back(0);
  // q.push_back(0);
  // q.push_back(-epsilon1);
////////////////////////////////////////////////////////////////////////////////////////



std::vector<long double> P_plus_LAB(const long double mjarg, const long double mkarg, const long double Mnarg, const std::vector<long double> &x){

  std::vector<long double> P;
  
  // 4 vector in S frame 
  std::vector<long double> P_S = P_plus_Sframe(mjarg, mkarg, Mnarg, x);

  // LAB transformation
  P = S_to_LAB(P_S, S_to_LAB_params(mjarg, mkarg, Mnarg, x) );

  // Check

  // if(Momentum_contains_nans(P_S)){

  //     std::cout<<"WARNING! PplusS momentum is nan-ed."<<std::endl; 
  //     std::cout<<"P: ("<<P_S[0]<<", "<<P_S[1]<<", "<<P_S[2]<<", "<<P_S[3]<<")"<<std::endl;
  // }

  // if(Momentum_contains_nans(P)){
  //     std::cout<<"WARNING! PplusLAB momentum is nan-ed."<<std::endl; 
  //     std::cout<<"P: ("<<P[0]<<", "<<P[1]<<", "<<P[2]<<", "<<P[3]<<")"<<std::endl;
  // }

  return P;
}

std::vector<long double> P_minus_LAB(const long double mjarg, const long double mkarg, const long double Mnarg, const std::vector<long double> &x){

  std::vector<long double> P;
  
  // 4 vector in S frame 
  std::vector<long double> P_S = P_minus_Sframe(mjarg, mkarg, Mnarg, x);

  // LAB transformation
  P = S_to_LAB(P_S, S_to_LAB_params(mjarg, mkarg, Mnarg, x) );

  // Check
  // if(Momentum_contains_nans(P_S)){

  //     std::cout<<"WARNING! PminusS momentum is nan-ed."<<std::endl; 
  //     std::cout<<"P: ("<<P_S[0]<<", "<<P_S[1]<<", "<<P_S[2]<<", "<<P_S[3]<<")"<<std::endl;
  // }

  // if(Momentum_contains_nans(P)){
  //     std::cout<<"WARNING! PminusLAB momentum is nan-ed."<<std::endl; 
  //     std::cout<<"P: ("<<P[0]<<", "<<P[1]<<", "<<P[2]<<", "<<P[3]<<")"<<std::endl;
  // }

  return P;
}