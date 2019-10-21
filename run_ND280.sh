#!/bin/sh

make

NEVENTS=100000
FLUXFILE=fluxes/ND280_horn250kA.dat
# RUN MU MU CHANNEL for all ND280 targets
./gen_SM -c 1 -z 6 -a 12 --fluxfile=$FLUXFILE --emin=0.01 --emax=15.0 -N $NEVENTS
./gen_SM -c 1 -z 8 -a 16 --fluxfile=$FLUXFILE --emin=0.01 --emax=15.0 -N $NEVENTS
./gen_SM -c 1 -z 82 -a 207 --fluxfile=$FLUXFILE --emin=0.01 --emax=15.0 -N $NEVENTS
./gen_SM -c 1 -z 1 -a 1 --pb --fluxfile=$FLUXFILE --emin=0.01 --emax=15.0 -N $NEVENTS
./gen_SM -c 1 -z 1 -a 1 --fluxfile=$FLUXFILE --emin=0.01 --emax=15.0 -N $NEVENTS