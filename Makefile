make: cross_sections.cxx integrator.cxx integrands.cxx observables.cxx integrator.cxx constants.h
	g++ gen_SM.cxx integrator.cxx integrands.cxx cross_sections.cxx observables.cxx constants.cxx -o gen_SM -std=c++11 -lm -lcuba -I ~/local/include/ -L ~/local/lib/ -lboost_program_options
