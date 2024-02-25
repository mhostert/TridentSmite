#include "observables.h"
  
/*

Dealing with 4 vectors in two frames: S and Lab frame.

S frame is defined as the center of mass frame of the neutrino and the photon,
More precisely, it is  k1_vec + q_vec = 0

*/

std::vector<long double> S_to_LAB(const std::vector<long double> &x, const std::vector<long double> &params){

  std::vector<long double> x_L(4);

  long double gamma = params[0];
  long double beta = params[1];
  long double sina = params[2];
  long double cosa = params[3];
  long double sind = params[4];
  long double cosd = params[5];

  x_L[0] = gamma*x[0] - gamma*beta*(sind*x[1] - cosd*x[3])    ;
  x_L[1] = gamma*beta*sina*x[0] + (cosa*cosd -gamma*sina*sind)*x[1] + (cosa*sind+sina*cosd*gamma)*x[3]  ;
  x_L[2] = x[2] ;
  x_L[3] = gamma*beta*cosa*x[0] - (sina*cosd + gamma*cosa*sind)*x[1] - (sina*sind-gamma*cosa*cosd)*x[3]  ;

  return x_L;
}

std::vector<long double> S_to_LAB_params(const long double mjarg, const long double mkarg, const long double Mnarg, const std::vector<long double> &x){

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
}

std::vector<std::vector<long double>> P_Sframe(const long double mjarg, const long double mkarg, const long double Mnarg, const std::vector<long double> &x){

  long double Q2   = x[0];      // -q^2 x1
  long double dot_q1   = x[1];      // q.k1 x2
  long double dot_qplus    = x[2]; // q.p_plus x3
  long double dot_qminus   = x[3]; // q.p_minus x4
  long double dot_1plus    = x[4]; // k1.p_plus x5
  long double dot_1minus   = x[5]; // k1.p_minus x6
  long double phi_minus = x[6]; // x7 random number (0,2pi) for the special-frame azimuthal angle (taken to be the l_minus azimuthal angle in that frame).
  // long double phi = x[7]; // x8 random number (0,2pi) for the trivial overall azimuthal angle
  // long double Enu   = x[8];     // Enu in lab

  /************************************************************************************************/
  long double m_minus = mjarg;
  long double m_plus  = mkarg;

  long double shat = 2*dot_q1 - Q2; // (k1 + q) COM energy
  long double dot_12 = dot_q1 - dot_1plus - dot_1minus; // k1.p2
  long double dot_q2 = dot_q1 - dot_qplus - dot_qminus - Q2; // q.p2
  long double Eplus = (dot_1plus + dot_qplus)/sqrt(shat);
  long double Eminus = (dot_1minus + dot_qminus)/sqrt(shat);
  long double Ep2 = (dot_12+dot_q2)/sqrt(shat);

  long double Pplus = sqrt(Eplus*Eplus - m_plus*m_plus);
  long double Pminus = sqrt(Eminus*Eminus - m_minus*m_minus);
  long double Pp2 = Ep2;

  long double Eq = (dot_q1-Q2)/(sqrt(shat)); // energy of q in S frame
  long double p1 = (dot_q1)/(sqrt(shat)); // momentum of k1 and -q in S frame

  long double cosplus = (p1*dot_qplus - Eq*dot_1plus)/(p1*Pplus*sqrt(shat));
  long double sinplus = sqrt(1-cosplus*cosplus);
  long double cosminus = (p1*dot_qminus - Eq*dot_1minus)/(p1*Pminus*sqrt(shat));
  long double sinminus = sqrt(1-cosminus*cosminus);
  
  long double cosp2 = (p1*dot_q2 - Eq*dot_12)/(p1*Pp2*sqrt(shat));
  long double sinp2 = sqrt(1-cosp2*cosp2);

  // We have to work a bit harder to get the azmithual angle of phi_plus from phi_minus.
  // The z_ij = dot_ij 
  long double dot_pm = (-1)*(dot_q1 - dot_qplus - dot_qminus - dot_1plus - dot_1minus + 0.5*(m_minus*m_minus + m_plus*m_plus - Q2));
  long double dot_m2 = (dot_q1 - dot_qplus - dot_1plus - 0.5*(m_minus*m_minus - m_plus*m_plus + Q2));
  /************************************************************************************************/

  // The sign of p2 does not matter it seems.
  long double SIGN_PM = +1;
  long double SIGN_M2 = -1;
  long double phi_plus = phi_minus + SIGN_PM * acos((Eplus*Eminus - dot_pm - Pplus*Pminus*cosplus*cosminus)/(Pplus*Pminus*sinplus*sinminus));
  long double phi_p2 = phi_minus + SIGN_M2 * acos((Ep2*Eminus - dot_m2 - Pp2*Pminus*cosp2*cosminus)/(Pp2*Pminus*sinp2*sinminus));
  
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
}

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

  // P1.push_back(p1);
  // P1.push_back(0.0);
  // P1.push_back(0.0);
  // P1.push_back(p1);

  // // NoP4 P4e need to boost it all back into the lab frame.
  // PN_L.push_back(mN);
  // PN_L.push_back(0.0);
  // PN_L.push_back(0.0);
  // PN_L.push_back(0.0);

  // q_L.push_back(q0);
  // q_L.push_back(qP*sinq);
  // q_L.push_back(0);
  // q_L.push_back(qP*cosq);

  // q.push_back(Eq);
  // q.push_back(0);
  // q.push_back(0);
  // q.push_back(-p1);
////////////////////////////////////////////////////////////////////////////////////////


  std::vector<std::vector<long double>> P_LAB(const long double mjarg, const long double mkarg, const long double Mnarg, const std::vector<long double> &x){

  std::vector<std::vector<long double>> P_LAB;
  
  // 4 vectors in S frame 
  std::vector<std::vector<long double>> vector_of_P_S = P_Sframe(mjarg, mkarg, Mnarg, x);

  // LAB transformation
  for (int i = 0; i < vector_of_P_S.size(); i++)
  {
    P_LAB.push_back( S_to_LAB(vector_of_P_S[i], S_to_LAB_params(mjarg, mkarg, Mnarg, x) ) );
  }

  return P_LAB;
}