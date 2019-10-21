#  SM/BSM Neutrino Trident Event Generator 


    #######################################################
    # ___ ____ _ ___  ____ _  _ ___ ____ _  _ _ ___ ____  #
    #  |  |__/ | |  \ |___ |\ |  |  [__  |\/| |  |  |___  #
    #  |  |  \ | |__/ |___ | \|  |  ___] |  | |  |  |___  #
    #                                                     #
    #                                                     # 
    #                        #########                    # 
    #                      ##                             # 
    #           ######################                    # 
    #                      ##                             #
    #                        #########                    #
    #                                                     #
    #######################################################


Generates neutrino trident events based on the calculation in the following papers:

https://arxiv.org/abs/1807.10973
https://arxiv.org/abs/1902.08579

Requirements
- working installation of the CUBA integration library (see http://www.feynarts.de/cuba/).
- C++ boost (-lboost_program_options)

See Makefile for environment variables.

## Usage

After your own setup of the code, you can run

```
	make 
	./gen_SM --help

	Allowed options:
	  -h [ --help ]                         produce help message
	  -N [ --nevents ] arg (=10000)         Number of HEPevt events to generate
	  -c [ --channel ] arg (=1)             Trident channel to use (see README for 
	                                        definition)
	  -z [ --znumber ] arg (=1)             Target proton number
	  -a [ --anumber ] arg (=1)             Target mass number (e.g. A = 12 for 
	                                        Carbon)
	  -b [ --pb ]                           Pauli blocking (only include if 
	                                        computing scattering on bound protons)
	  -l [ --emin ] arg (=0.100000000000000005551)
	                                        Minimum Enu to sample from
	  -u [ --emax ] arg (=40)               Maximum Enu to sample from
	  -f [ --fluxfile ] arg (=fluxes/uniform_0.1_200_GeV.dat)
	                                        Neutrino flux file to use
	  -m [ --mzprime ] arg (=1)             Zprime mass
	  -g [ --gprime ] arg (=0)              Zprime coupling

```

A few neutrino fluxes are available under "fluxes/", including files used for a uniform neutrino energy distribution. By specifying 'emin' and 'emax', one can sample only that energy window of the flux.

Whenever A>1, a coherent cross section is calculated. For A=1, if Z=0 (Z=1), then elastic scattering on neutrons (protons) is calculated. If option -b(--pb) is set Pauli blocking effects are taken into account, otherwise these are ignored.

CHANNEL stands for the desired trident channel where for the channel nu_a -> nu_b ell_c^+ ell_d^- = abcd, we have:

```	
	eeee 0
	mmmm 1
	emme 2
	mmee 3
	eemm 4
	meem 5
	eett 6
	mmtt 7
	tttt 8
	ette 9
	mttm 10
	teet 11
	tmmt 12
	ttee 13
	ttmm 14
```