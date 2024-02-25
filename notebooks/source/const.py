import numpy as np
import scipy
from scipy import interpolate

## MASSES in GeV
rad_to_deg = 180.0/np.pi

m_proton = 0.93827208816 # GeV
m_neutron = 0.93956542052 # GeV
m_AVG = (m_proton+m_neutron)/2. # GeV

mW = 80.37912 # GeV
mZ = 91.187621 # GeV

GeV2_to_cm2 = 3.89379372e-28 # hbar c = 197.3269804e-16 GeV.cm

Gf=1.16637876e-5 # Fermi constant (GeV^-2)

alphaQED = 1.0/137.03599908421 # Fine structure constant at q2 -> 0

m_e =  0.5109989500015e-3 # GeV
m_mu =  0.1134289257 # GeV
m_tau =  1.77682 # GeV

#Avogadro's number
NAvo = 6.02214076*1e23

# speed of light (PDG) m/s
c_LIGHT = 299792458


# constants for normalization
invm2_to_incm2=1e-4
zb_to_cm2 = 1e-45
nucleons_to_tons = NAvo*1e6/m_AVG


# trident channel definitions
eeee = 0
mmmm = 1
emme = 2
mmee = 3
eemm = 4
meem = 5
eett = 6
mmtt = 7
tttt = 8
ette = 9
mttm = 10
teet = 11
tmmt = 12
ttee = 13
ttmm = 14

