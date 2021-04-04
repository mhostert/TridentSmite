#!/bin/sh

make

NEVENTS=100000
FLUXFILE=fluxes/SBN/Flux_SBN_MB.dat
EMAX=5.0

# RUN MU MU CHANNEL for all DUNE ND targets
./gen_SM -c 1 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=0.01 --emax=$EMAX -N $NEVENTS
./gen_SM -c -1 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=0.01 --emax=$EMAX -N $NEVENTS

./gen_SM -c 4 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=0.01 --emax=$EMAX -N $NEVENTS
./gen_SM -c -4 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=0.01 --emax=$EMAX -N $NEVENTS

./gen_SM -c 1 -z 1 -a 1 --pb --fluxfile=$FLUXFILE --emin=0.01 --emax=$EMAX -N $NEVENTS
./gen_SM -c -1 -z 1 -a 1 --pb --fluxfile=$FLUXFILE --emin=0.01 --emax=$EMAX -N $NEVENTS
