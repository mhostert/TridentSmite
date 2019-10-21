#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
/* ********************************************************* */
// Physical variables

#define eeee_flag 0
#define mmmm_flag 1
#define emme_flag 2
#define mmee_flag 3
#define eemm_flag 4
#define meem_flag 5

#define eett_flag 6
#define mmtt_flag 7
#define tttt_flag 8
#define ette_flag 9
#define mttm_flag 10
#define teet_flag 11
#define tmmt_flag 12
#define ttee_flag 13
#define ttmm_flag 14

// flags for the object scan
#define WFLUX 0
#define WOFLUX 1

#define NO_SAMPLES 0
#define PRINT_SAMPLES 1

#define SMonly 0
#define INTERFERENCE 1
#define BSMonly 2
#define SMandBSM 3

#define NO_BLOCKING 0
#define W_BLOCKING 1


// Flags
#define e_flag 0
#define mu_flag 1
#define tau_flag 2

// BSM Flags
#define ZPRIME 0
#define WPRIME 1
#define Q2CUT 2

// Masses
#define m_e 0.5110000000000000000000e-3 // GeV
#define m_mu 105.70000000000000000000e-3 // GeV
#define m_tau 1.7770000000000000000000 // GeV

#define m_pi0 0.13497 // GeV
#define m_pic 0.13957 // GeV

#define m_proton 0.93827 // GeV
#define m_neutron 0.93956 // GeV
#define m_AVG 0.9389 // GeV

#define mW 80.385 // GeV
#define mZ 91.187 // GeV

// Useful constants
#define Gf 1.1663787e-5 // Fermi constant (GeV^-2)
#define alpha_QED 1.0/137.035 // Fine structure constant
#define alphaQED 1.0/137.035 // Fine structure constant

#define sw2 0.2223 // sin of Weinber angle squared
#define gweak 8*Gf/sqrt(2)*mW*mW
#define gL sw2 - 1.0/2.0
#define gR sw2
#define Ca gR - gL
#define Cv gR + gL + 1.0

#define Vud 0.97425 // CKM matrix first entry |Vud| w mod

#define f_decay_pi0 0.130 // Neutral pion decay constant (GeV)
#define f_decay_pic  0.13041 // Charged pion decay constant (GeV)


// FORM FACTOR DEFINITIONS

#define LAMBDA_QCD 0.217 // GeV
#define Q_MAX LAMBDA_QCD/pow(A,1.0/3.0)// + 0.2*LAMBDA_QCD/pow(A,1.0/3.0) //GeV
#define Q_MAX_DIFF 1.0// + 0.2 //GeV


#define M_avg (m_proton+m_neutron)/2.0
#define qsi 4.7 // (mu_p - mu_m)/mu_N 
#define MV 0.84 // GeV (Vectorial mass)
#define MAG_MOM_N -1.913 // in units of the nuclear magneton e hbar / 2 mp  
#define MAG_MOM_P 2.792 // in units of the nuclear magneton e hbar / 2 mp 

#define P_FERMI 0.235 // GeV fermi momentum of nucleus


#define N_discrete 100
#define MAX_ITER_LIM 1e7

#define GeV2_to_cm2 3.9204e-28 // 1 GeV ^-2 = .... cm^2  

inline long double SQR(long double x) { return x*x; }

inline long double heaviside(long double x) { return 1.0*(x > 0.0); }

									// x = |Q|
inline long double G_dip(long double x) { return SQR( 1.0 / (1.0 + SQR(x/MV) ) ); }// Recall q = sqrt(-q^2) = |Q^2|, so tau = -|Q|^2/4M^2

/* ********************************************************* */


/* ********************************************************* */
// Tedious variables for the integration algorithm
#define NDIM        9
#define NCOMP 		1
#define USERDATA 	NULL
#define NVEC 		1  
#define EPSREL 		5e-4
#define EPSABS 		5e-4

#define VERBOSE 	0
#define LAST 		1
#define SEED 	  time(0)
#define MINEVAL 	4e5
#define MAXEVAL     4e5
#define NSTART 		1e3
#define NINCREASE 	1e3
#define NBATCH 		1
#define GRIDNO 		0
#define STATEFILE 	NULL
#define SPIN 		NULL

#define SKIP_PRINTS 0.9 // print integrand after SKIP_PRINTS of evaluations

// SUAVE ARGUMENTS --> Not tested yet (14/12)
#define NNEW 1e4
#define NMIN 2
#define FLATNESS 20

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

#define D_KEY1 47
#define D_KEY2 40
#define D_KEY3 40
#define D_MAXPASS 5
#define D_BORDER 0.
#define D_MAXCHISQ 100.0
#define D_MINDEVIATION .25
#define D_REGIONS 50


/* ********************************************************* */
// RANDOM NUMBER GENERATORS
/////////////////////////////////////////////////////
long double UniformRand();
long double NormalRand(long double mean, long double stddev);

// Does it have nans?
int Momentum_contains_nans(const std::vector<long double> &P);

// strings and printing
void pretty_print(std::string CHANNEL);
std::vector<std::string> list_of_trident_channels(std::vector<std::string> &sList);
std::string getLastLine(std::ifstream& in);


class trident_channel{
  // What channel is being simulated nu_alpha  -> nu l2^+ l1^- 
 	public:
 	int nu_alpha,l1,l2;
 	std::string channel_name;
 	trident_channel(int C); 
};

/* ********************************************************* */
#endif