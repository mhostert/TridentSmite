import numpy as np
from . import const
from . import gen
from copy import copy

def set_new_physics_for_all(gprimeV,gprimeA,mzprime):
	MINOS_LE_FHC.set_new_physics(gprimeV,gprimeA,mzprime)
	MINOS_LE_RHC.set_new_physics(gprimeV,gprimeA,mzprime)
	MINOS_ME_FHC.set_new_physics(gprimeV,gprimeA,mzprime)
	NOVA_RHC.set_new_physics(gprimeV,gprimeA,mzprime)
	NOVA_FHC.set_new_physics(gprimeV,gprimeA,mzprime)
	ND280_RHC.set_new_physics(gprimeV,gprimeA,mzprime)
	ND280_FHC.set_new_physics(gprimeV,gprimeA,mzprime)


#######################
# MINOS
#######################

## relative mass composition given in https://arxiv.org/pdf/1902.00558.pdf
## Low energy NuMI +  FHC mode
MINOS_LE_FHC = gen.exp_setup("MINOS_LE_FHC")
MINOS_LE_FHC.fluxfile="fluxes/NUMI_FHC_LE.dat"
MINOS_LE_FHC.channels=np.array([const.mmmm, -const.mmmm])
MINOS_LE_FHC.Erange = np.array([0.01,50.0])
MINOS_LE_FHC.totmass= 28.6
MINOS_LE_FHC.Z =    np.array([6,26])
MINOS_LE_FHC.A =    np.array([12,56])
MINOS_LE_FHC.mass = np.array([0.2*MINOS_LE_FHC.totmass, 0.8*MINOS_LE_FHC.totmass])
MINOS_LE_FHC.POT  =  10.56e20
MINOS_LE_FHC.norm = MINOS_LE_FHC.POT*const.nucleons_to_tons*const.zb_to_cm2*const.invm2_to_incm2/1.0


## Low energy NuMI +  FHC mode
MINOS_LE_RHC = gen.exp_setup("MINOS_LE_RHC")
MINOS_LE_RHC.fluxfile="fluxes/NUMI_RHC_LE.dat"
MINOS_LE_RHC.channels=np.array([const.mmmm, -const.mmmm])
MINOS_LE_RHC.Erange = np.array([0.01,50.0])
MINOS_LE_RHC.totmass= 28.6
MINOS_LE_RHC.Z =    np.array([6,26])
MINOS_LE_RHC.A =    np.array([12,56])
MINOS_LE_RHC.mass = np.array([0.2*MINOS_LE_RHC.totmass, 0.8*MINOS_LE_RHC.totmass])
MINOS_LE_RHC.POT  =  3.36e20
MINOS_LE_RHC.norm = MINOS_LE_RHC.POT*const.nucleons_to_tons*const.zb_to_cm2*const.invm2_to_incm2/1.0

## Low energy NuMI +  FHC mode
MINOS_ME_FHC = gen.exp_setup("MINOS_ME_FHC")
MINOS_ME_FHC.fluxfile="fluxes/NUMI_FHC_ME_unofficial.dat"
MINOS_ME_FHC.channels=np.array([const.mmmm])
MINOS_ME_FHC.Erange = np.array([0.01,15.0])
MINOS_ME_FHC.totmass= 28.6
MINOS_ME_FHC.Z =    np.array([6,26])
MINOS_ME_FHC.A =    np.array([12,56])
MINOS_ME_FHC.mass = np.array([0.2*MINOS_ME_FHC.totmass, 0.8*MINOS_ME_FHC.totmass])
MINOS_ME_FHC.POT  =  9.69e20
MINOS_ME_FHC.norm = MINOS_ME_FHC.POT*const.nucleons_to_tons*const.zb_to_cm2*const.invm2_to_incm2/1.0


#######################
# NOvA
#######################
# masses in t 
totmass=193
massTi=3.2e-2*totmass
massCl=16.1e-2*totmass
massO=3.00e-2*totmass
massC=66.7e-2*totmass
massH=10.8e-2*totmass

NOVA_FHC = gen.exp_setup("NOVA_FHC")
NOVA_FHC.fluxfile="fluxes/NOvA_FHC.dat"
NOVA_FHC.channels=np.array([const.mmmm,-const.mmmm])
NOVA_FHC.Erange = np.array([0.01,20.0])
NOVA_FHC.totmass= totmass # t
NOVA_FHC.Z =    np.array([6,17,22,8,1])
NOVA_FHC.A =    np.array([12,35,48,16,1])
NOVA_FHC.mass = np.array([massC,massCl,massTi,massO,massH])
NOVA_FHC.POT  = 1.36e15 #1e6
NOVA_FHC.norm = const.nucleons_to_tons*const.zb_to_cm2*NOVA_FHC.POT*const.invm2_to_incm2/1.0

NOVA_RHC = gen.exp_setup("NOVA_RHC")
NOVA_RHC.fluxfile="fluxes/NOVA_RHC.dat"
NOVA_RHC.channels=np.array([const.mmmm,-const.mmmm])
NOVA_RHC.Erange = np.array([0.01,20.0])
NOVA_RHC.totmass= totmass # t
NOVA_RHC.Z =    np.array([6,17,22,8,1])
NOVA_RHC.A =    np.array([12,35,48,16,1])
NOVA_RHC.mass = np.array([massC,massCl,massTi,massO,massH])
NOVA_RHC.POT  = 1.25e15 #1e6
NOVA_RHC.norm = const.nucleons_to_tons*const.zb_to_cm2*NOVA_RHC.POT*const.invm2_to_incm2/1.0


#######################
# ND280 T2K
#######################
# masses in t 
massZn=0.8
massCu=0.4
massO=3.32
massC=5.4
massPb=12.0
massH=0.42

dE=0.05

ND280_FHC = gen.exp_setup("ND280_FHC")
ND280_FHC.fluxfile="fluxes/ND280_FHC.dat"
ND280_FHC.channels=np.array([const.mmmm,-const.mmmm])
ND280_FHC.Erange = np.array([0.01,20.0])
ND280_FHC.Z =    np.array([82,29,30,8,6,1])
ND280_FHC.A =    np.array([207,65,65,16,12,1])
ND280_FHC.mass = np.array([massPb,massCu,massZn,massO,massC,massH])
ND280_FHC.totmass= np.sum(ND280_FHC.mass)
ND280_FHC.POT  = 1.97 #1e21
ND280_FHC.norm = const.nucleons_to_tons*const.zb_to_cm2*ND280_FHC.POT/dE

ND280_RHC = gen.exp_setup("ND280_RHC")
ND280_RHC.fluxfile="fluxes/ND280_RHC.dat"
ND280_RHC.channels=np.array([const.mmmm,-const.mmmm])
ND280_RHC.Erange = np.array([0.01,20.0])
ND280_RHC.Z =    np.array([82,29,30,8,6,1])
ND280_RHC.A =    np.array([207,65,65,16,12,1])
ND280_RHC.mass = np.array([massPb,massCu,massZn,massO,massC,massH])
ND280_RHC.totmass= np.sum(ND280_RHC.mass)
ND280_RHC.POT  = 1.63 #1e21
ND280_RHC.norm = const.nucleons_to_tons*const.zb_to_cm2*ND280_RHC.POT/dE


#######################
# CHARM-II
#######################
# masses in t 
massZn=0.8

dE=1.0

CHARMII_FHC = gen.exp_setup("CHARMII_FHC")
CHARMII_FHC.fluxfile="fluxes/CHARMII_FHC.dat"
CHARMII_FHC.channels=np.array([const.mmmm,-const.mmmm])
CHARMII_FHC.Erange = np.array([1.5,200])
CHARMII_FHC.Z =    np.array([11])
CHARMII_FHC.A =    np.array([20.7])
CHARMII_FHC.mass = np.array([574])
CHARMII_FHC.totmass= np.sum(CHARMII_FHC.mass)
CHARMII_FHC.POT  = 1.5e19 # POTs
CHARMII_FHC.norm = const.nucleons_to_tons*const.zb_to_cm2*CHARMII_FHC.POT/dE
