#!/bin/sh

cd ..

NEVENTS=100000
FLUXFILE=fluxes/SBN/Flux_SBN_MB.dat

# create a folder with the name of this current file to store events
THIS_FILE_NO_EXT="$(basename -- $0)"
EVENTS_PATH="data/${THIS_FILE_NO_EXT%.*}/"
echo $EVENTS_PATH
EMAX=5.00
EMIN=0.001

# RUN MU MU CHANNEL for all ND280 targets
./gen_SM -c 1 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2
./gen_SM -c -1 -z 18 -a 40 --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2

# p elastic
./gen_SM -c 1 -z 1 -a 1 --pb --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2
./gen_SM -c -1 -z 1 -a 1 --pb --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2

# n elastic
./gen_SM -c 1 -z 0 -a 1 --pb --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2
./gen_SM -c -1 -z 0 -a 1 --pb --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2