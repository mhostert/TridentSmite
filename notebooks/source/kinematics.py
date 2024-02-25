import numpy as np
import scipy
import re
import os

from . import Cfv
from . import const

# Class with all the relevant quantities of a given trident channel/mode
class kin:
	##############################################    
	# Import datafile and set lepton masses
	def __init__(self, dic_of_events, mplus, mminus, Emin=0.0, Emax=20.0):
		
		# events
		self.dic = dic_of_events
		
		# charge lepton masses 
		self.mp = mplus
		self.mm = mminus

		# charge lepton masses 
		self.Enu_min = Emin
		self.Enu_max = Emax

		# event properties
		self.Z = self.dic['Z'] 
		self.A = self.dic['A'] 
		self.PB = self.dic['PB'] 
		self.cp = self.dic['cp'] 

		# relevant quantities -- unfortunally my Cfv script only takes doubles not long doubles
		self.Enu = self.dic['Enu'].astype(np.float64) 
		self.Q2  = self.dic['Q2'].astype(np.float64) 
		self.w   = self.dic['weights'].astype(np.float64)
		self.Pnu = self.dic['Pnu'].astype(np.float64).T
		self.Pplus = self.dic['Pplus'].astype(np.float64).T
		self.Pminus = self.dic['Pminus'].astype(np.float64).T

		print("Total samples imported: ", np.size(self.w))

		# mask of all events

		self.all_events = np.ones(np.size(self.w)).astype(bool)
		
		self.only_coh = (self.PB == False)
		self.only_incoh = (self.PB == True)

		self.only_nu = (self.cp == 1)
		self.only_nubar = (self.cp == -1)


	##############################################
	# Function to reweight the neutrino flux
	# the flux files has to have the usual format of 7 columns:
	#
	# Enu(GeV) nue numu nutau nue_bar numu_bar nutau_bar
	#
	# The uniform flux is used in the MC is set to 1, so user can choose the units of the new weights
	def Reweight(self, old_flux_file, new_flux_file, flavour):
	    M_old = np.loadtxt(old_flux_file, dtype=np.float64)
	    en_old = M_old[:,0]
	    f_old  = M_old[:,flavour]     
	    M = np.loadtxt(new_flux_file, dtype=np.float64)
	    en = M[:,0]
	    f  = M[:,flavour]
	    # interpolate numu flux and fill out of range values w/ 0
	    flux_old = interp1d(en_old, f_old, fill_value=0, bounds_error=False)
	    flux = interp1d(en, f, fill_value=0, bounds_error=False)
	    
	    self.rw = self.w * np.where( flux(self.Enu)*flux_old(self.Enu) > 0, flux(self.Enu)/flux_old(self.Enu) , 0)

	##############################################    
	# flux convolved cross section calculated from the samples in 1e-48cm^2*arbitrary_flux_units (using uniform flux)
	def total_fxs(self):
	    return np.sum(self.w)

	# flux convolved cross section calculated from the samples in 1e-48cm^2*arbitrary_flux_units (using new flux)
	def total_fxs_new_flux(self):
	    return np.sum(self.rw)
	##############################################    


	##############################################    
	# invariant mass of the two charged leptons 
	def invmassSQR(self, mask=...):
	    return Cfv.dot4(self.Pplus[mask]+self.Pminus[mask], self.Pplus[mask]+self.Pminus[mask])

	def invmass(self, mask=...):
	    return np.sqrt(Cfv.dot4(self.Pplus[mask]+self.Pminus[mask], self.Pplus[mask]+self.Pminus[mask]))

	# separation angle of the two charged leptons in rads 
	def sepangle(self, mask=...):
	    return np.arccos(Cfv.dot3(self.Pplus[mask], self.Pminus[mask])/np.sqrt(Cfv.dot3(self.Pplus[mask], self.Pplus[mask])*Cfv.dot3(self.Pminus[mask], self.Pminus[mask])))

	# angle of wrt to the beam of negatively charged lepton in rads 
	def minus_beam(self, mask=...):
	    return np.arccos(self.Pminus[:,3][mask]/np.sqrt(Cfv.dot3(self.Pminus[mask], self.Pminus[mask])))

	# angle of wrt to the beam of positively charged lepton in rads 
	def plus_beam(self, mask=...):
	    return np.arccos(self.Pplus[:,3][mask]/np.sqrt(Cfv.dot3(self.Pplus[mask], self.Pplus[mask])))

	# Energy of negatively charged lepton in GeV
	def Eminus(self, mask=...):
	    return self.Pminus[:,0][mask]

	# Energy of positively charged lepton in GeV
	def Eplus(self, mask=...):
	    return self.Pplus[:,0][mask]

	# angle_of_cone_wrt_beam 
	def angle_of_cone_wrt_beam(self, mask=None):
	    return np.arccos((self.Pplus[mask]+self.Pminus[mask])[:,3]/np.sqrt(Cfv.dot3(self.Pplus[mask]+self.Pplus[mask],self.Pplus[mask]+self.Pplus[mask])))

