#!/bin/sh

make

NEVENTS=100000
FLUXFILE=fluxes/DUNE_ND_neutrino_Nov2017.dat

# RUN MU MU CHANNEL for all DUNE ND targets
./gen_SM -c 1 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=0.01 --emax=30.0 -N $NEVENTS
# ./gen_SM -c -1 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=0.01 --emax=30.0 -N $NEVENTS
# ./gen_SM -c 4 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=0.01 --emax=30.0 -N $NEVENTS
# ./gen_SM -c -4 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=0.01 --emax=30.0 -N $NEVENTS
# ./gen_SM -c 1 -z 1 -a 1 --pb --fluxfile=$FLUXFILE --emin=0.01 --emax=30.0 -N $NEVENTS
# ./gen_SM -c -1 -z 1 -a 1 --pb --fluxfile=$FLUXFILE --emin=0.01 --emax=30.0 -N $NEVENTS
# ./gen_SM -c 1 -z 0 -a 1 --pb --fluxfile=$FLUXFILE --emin=0.01 --emax=30.0 -N $NEVENTS


# ./gen_SM -c 1 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=0.01 --emax=30.0 -N $NEVENTS --gprimeV=0.00098 --mzprime=0.1
# ./gen_SM -c 1 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=0.01 --emax=30.0 -N $NEVENTS --gprimeV=0.00056 --mzprime=0.015
# ./gen_SM -c 1 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=0.01 --emax=30.0 -N $NEVENTS --gprimeV=0.1 --mzprime=0.015