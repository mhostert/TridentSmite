#ifndef FLUX_H_
#define FLUX_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>



#define DET_NOCUTS 0
#define DET_SBND 1
#define DET_MUBOONE 2
#define DET_ICARUS 3
#define DET_PS191 4
#define DET_ALL_SBN 5
#define DET_DUNE_RHC_ND 6
#define DET_DUNE_FHC_ND 7
#define DET_NUSTORM_LAr_ND 8
#define DET_T2K_INGRID 9
#define DET_SHiP 10

#define MPION  0.13957
#define MPI0   0.13497
#define MKAON  0.49367
#define MMU    0.10566
#define	ME     0.00051
#define M2GEV 5.06842e15
#define INVS2GEV 6.58e-25
#define FKA 0.1561
#define VUS 0.2252

#define FUDGE_SCALE  false
#define FUDGE_SHIFT false
;

double heaviside(double);
double geometry_factor(double);


double fluxDUNE_FD_e(double E);
double fluxDUNE_FD_mu(double E);

double fluxDUNE_ND_e(double E);
double fluxDUNE_ND_mu(double E);
double fluxDUNE_ND_ebar(double E);
double fluxDUNE_ND_mubar(double E);

double fluxDUNE_RHC_ND_e(double E);
double fluxDUNE_RHC_ND_mu(double E);
double fluxDUNE_RHC_ND_ebar(double E);
double fluxDUNE_RHC_ND_mubar(double E);

double fluxSBND_e(double E);
double fluxSBND_mu(double E);

double fluxMB_e(double E);
double fluxMB_mu(double E);

double fluxT2K_nue(double E);
double fluxT2K_nuebar(double E);
double fluxT2K_numu(double E);
double fluxT2K_numubar(double E);

double fluxNuSTORM_mudecay_numubar(double E);
double fluxNuSTORM_mudecay_nue(double E);
double fluxNuSTORM_pidecay_numu(double E);

double fluxSHiP_nue_nuebar(double E);
double fluxSHiP_numu_numubar(double E);
double fluxSHiP_nutau_nutaubar(double E);


#endif

