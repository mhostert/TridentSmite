#!/bin/sh

# make
cd ..
NEVENTS=100000
FLUXFILE=fluxes/DUNE_ND_neutrino_Nov2017.dat

EMIN=0.1
EMAX=30.0


# Argon40
./gen_SM -c 1 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --gprimeV=$1 --mzprime=$2 --text
./gen_SM -c -1 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --gprimeV=$1 --mzprime=$2 --text

# proton
./gen_SM -c 1 -z 1 -a 1 --pb --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --gprimeV=$1 --mzprime=$2 --text
./gen_SM -c -1 -z 1 -a 1 --pb --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --gprimeV=$1 --mzprime=$2 --text

# neutron
./gen_SM -c 1 -z 0 -a 1 --pb --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --gprimeV=$1 --mzprime=$2 --text
./gen_SM -c -1 -z 0 -a 1 --pb --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --gprimeV=$1 --mzprime=$2 --text