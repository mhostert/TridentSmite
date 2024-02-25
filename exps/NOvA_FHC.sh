#!/bin/sh

cd ..

NEVENTS=100000
FLUXFILE=fluxes/NOvA_FHC.dat

# create a folder with the name of this current file to store events
THIS_FILE_NO_EXT="$(basename -- $0)"
EVENTS_PATH="data/${THIS_FILE_NO_EXT%.*}/"
echo $EVENTS_PATH
EMAX=20.00
EMIN=0.01

# coh
./gen_SM -c 1 -z 6 -a 12 --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2
./gen_SM -c -1 -z 6 -a 12 --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2

./gen_SM -c 1 -z 8 -a 16 --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2
./gen_SM -c -1 -z 8 -a 16 --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2

./gen_SM -c 1 -z 17 -a 35 --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2
./gen_SM -c -1 -z 17 -a 35 --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2

./gen_SM -c 1 -z 22 -a 48 --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2
./gen_SM -c -1 -z 22 -a 48 --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2

# p elastic
./gen_SM -c 1 -z 1 -a 1 --pb --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2
./gen_SM -c -1 -z 1 -a 1 --pb --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2

# n elastic
./gen_SM -c 1 -z 0 -a 1 --pb --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2
./gen_SM -c -1 -z 0 -a 1 --pb --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2

# hydrogen
./gen_SM -c 1 -z 1 -a 1 --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2
./gen_SM -c -1 -z 1 -a 1 --fluxfile=$FLUXFILE --emin=$EMIN --emax=$EMAX -N $NEVENTS --path $EVENTS_PATH --gprimeV=$1 --mzprime=$2