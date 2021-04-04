import numpy as np
import scipy
from scipy import interpolate

## FLAGS
electron = 1
muon = 2
tau = 3

COHLH=1
DIFLH=2
COHRH=3
DIFRH=4

DIRAC = 'dirac'
MAJORANA = 'majorana'

HM = 'heavy_mediator'
LM = 'light_mediator'


# fix experiment to use
THREEPLUSONE = 89
THREEPLUSTWO = 90

charged_pion = 400
charged_kaon = 500

neutral_pion = 411
neutral_kaon = 511
neutral_eta = 611

neutrino0 = 99
neutrino1 = 11
neutrino2 = 22
neutrino3 = 33
neutrino4 = 44
neutrino5 = 55

neutrino_light 	= 999
neutrino_electron 	= 10
neutrino_muon 		= 20
neutrino_tau 		= 30
neutrino_sterile	= 40
neutrino_dark 		= 50

## MASSES in GeV

higgsvev = 246 # GeV

Me  = 511e-6 
Mmu = 0.105
Mtau = 1.777 

mproton = 0.938
mneutron = 0.939
MAVG = (mproton + mneutron)/2.0


Mw = 80.35 
Mz = 91
higgsvev = 246.22 

## COUPLINGS 
s2w = 0.231
sw = np.sqrt(0.231)
cw = np.sqrt(1.0 - s2w)

gl_lepton = -0.5 + s2w
gr_lepton = s2w

alphaQED = 1./137.0359991
eQED = np.sqrt(4.0*np.pi*alphaQED)

gvP = 1.0
Gf = 1.16e-5 # GeV^-2
g = np.sqrt(Gf*8/np.sqrt(2)*Mw*Mw)
gweak = g

# FORM FACTOR CONSTANTS
gA = 1.26
tau3 = 1

MAG_N = -1.913
MAG_P = 2.792

# charged hadrons
Mcharged_pion = 0.1396
Mcharged_kaon = 0.4937
Mcharged_rho = 0.7758

Fcharged_pion = 0.1307
Fcharged_kaon = 0.1598
Fcharged_rho = 0.220

# neutral hadrons
# charged hadrons
Mneutral_pion = 0.135
Mneutral_eta = 0.5478
Mneutral_rho = 0.7755

Fneutral_pion = 0.130
Fneutral_kaon = 0.1647
Fneutral_eta = 0.210



#########
# CKM elements
Vud = 0.97420
Vus = 0.2243
Vcd = 0.218
Vcs = 0.997
Vcb = 42.2e-3
Vub = 3.94e-3
Vtd = 8.1e-3 
Vts = 39.4e-3
Vtb = 1

################

#Avogadro's number
NAvo = 6.022*1e23
# from GeV^-2 to cm^2
GeV2_to_cm2 = 3.9204e-28

# speed of light (PDG) m/s
c_LIGHT = 299792458

## FORM FACTORS
def D(Q2):
	return 1.0/((1+Q2/0.84/0.84)**2)


def H1_p(Q2):
	tau = -Q2/4.0/mproton/mproton
	F1 = (D(Q2) - tau*MAG_P*D(Q2))/(1-tau)
	F2 = (MAG_P*D(Q2) - D(Q2))/(1-tau)
	return F1**2 - tau * F2**2
def H2_p(Q2):
	tau = -Q2/4.0/mproton/mproton
	F1 = (D(Q2) - tau*MAG_P*D(Q2))/(1-tau)
	F2 = (MAG_P*D(Q2) - D(Q2))/(1-tau)
	return (F1+F2)*(F1+F2)



def H1_n(Q2):
	tau = -Q2/4.0/mproton/mproton
	F1 = (- tau*MAG_N*D(Q2))/(1-tau)
	F2 = (MAG_N*D(Q2))/(1-tau)
	return F1**2 - tau * F2**2

def H2_n(Q2):
	tau = -Q2/4.0/mproton/mproton
	F1 = (-tau*MAG_N*D(Q2))/(1-tau)
	F2 = (MAG_N*D(Q2))/(1-tau)
	return (F1+F2)*(F1+F2)



def F1_EM(Q2):
	tau = -Q2/4.0/mproton/mproton
	return (D(Q2) - tau*MAG_P*D(Q2))/(1-tau)

def F2_EM(Q2):
	tau = -Q2/4.0/mproton/mproton
	return (MAG_P*D(Q2) - D(Q2))/(1-tau)



def F1_WEAK(Q2):
	tau = -Q2/4.0/mproton/mproton
	f = (0.5 - s2w)*(tau3)*(1-tau*(1+MAG_P-MAG_N))/(1-tau) - s2w*(1-tau*(1+MAG_P+MAG_N))/(1-tau)  
	return f*D(Q2)

def F2_WEAK(Q2):
	tau = -Q2/4.0/mproton/mproton
	f = (0.5 - s2w)*(tau3)*(MAG_P-MAG_N)/(1-tau) - s2w*(MAG_P+MAG_N)/(1-tau)  
	return f*D(Q2)

def F3_WEAK(Q2):
	f = gA*tau3/2.0/(1+Q2/1.02/1.02)**2
	return f


a = 0.523/0.197 # GeV^-1
def FEMcoh(Q,MA):
	r0 = 1.03*(MA**(1.0/3.0))/0.197 # GeV^-1
	return 3.0*np.pi*a/(r0**2 + np.pi**2 * a**2) * (np.pi*a *(1.0/np.tanh(np.pi*a*Q))*np.sin(Q*r0) - r0*np.cos(Q*r0))/(Q*r0*np.sinh(np.pi*Q*a))

def j1(z):
	return np.sin(z)/z/z - np.cos(z)/z

fm_to_GeV = 1.0/0.1975

def FWEAKcoh(Q,MA):
	a = 0.523/0.197 # GeV^-1
	s = 0.9*fm_to_GeV # fm to GeV^-1
	R = 3.9*fm_to_GeV*(MA/MAVG/40.0)**(1.0/3.0) # fm to GeV^-1
	return (3*np.abs(j1(Q*R)/Q/R))*np.exp(-Q*Q*s*s/2)

def Fpauli_blocking(Q2,MA):
	MA=12.0
	qvec = np.sqrt(Q2*(1+Q2/4.0/MA/MA))
	kf = 0.300
	return np.piecewise(qvec, [qvec <= 2*kf, qvec > 2*kf],[lambda x: 3.0/2.0*x/2.0/kf - 0.5*(x/2.0/kf)*(x/2.0/kf)*(x/2.0/kf), lambda x: 1])
