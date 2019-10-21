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


Generates trident events both in the Standard Model and with a Zprime.

Requires a working installation of the CUBA integration library (see http://www.feynarts.de/cuba/).

## Usage

After your own setup, you can run
	make 
	./gen_SM (int)CHANNEL (int)Znumber (int)Anumber (0 or 1)Block (long double)Enumin (long double)Enumax

where CHANNEL stands for the desired trident channel where for the channel nu_a -> nu_b ell_c^+ ell_d^- = abcd, we have:
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