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

Requires a working installation of the CUBA integration library (see http://www.feynarts.de/cuba/).

## Usage

After your own setup, you can run

```
	make 
	./gen_SM (int)CHANNEL (int)Znumber (int)Anumber (0 or 1)Block (long double)Enumin (long double)Enumax
```

A few neutrino fluxes are available under "fluxes/", including files used for a uniform neutrino energy distribution. By specifying 'Enumin' and 'Enumax', one can sample only that energy window of the flux.

Whenever A>1, a coherent cross section is calculated. For A=1, if Z=0 (Z=1), then elastic scattering on neutrons (protons) is calculated. If Block=1, Pauli blocking effects are taken into account, otherwise these are ignored.

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