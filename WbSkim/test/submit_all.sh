#!/bin/sh

VERSION=v01
CUT=0

if [ ! -z "$1" ]; then
  VERSION=$1
fi

if [ ! -z "$2" ]; then
  CUT=$2
fi

./submit.sh SingleElectron_2012A_22Jan13 $VERSION $CUT
./submit.sh SingleElectron_2012B_22Jan13 $VERSION $CUT
./submit.sh SingleElectron_2012C_22Jan13 $VERSION $CUT
./submit.sh SingleElectron_2012D_22Jan13 $VERSION $CUT

./submit.sh SingleMu_2012A_22Jan13 $VERSION $CUT
./submit.sh SingleMu_2012B_22Jan13 $VERSION $CUT
./submit.sh SingleMu_2012C_22Jan13 $VERSION $CUT
./submit.sh SingleMu_2012D_22Jan13 $VERSION $CUT

./submit.sh DYJetsToLL $VERSION $CUT
./submit.sh TTbar $VERSION $CUT
./submit.sh TToLeptons $VERSION $CUT
./submit.sh TBarToLeptons $VERSION $CUT
./submit.sh WW $VERSION $CUT
./submit.sh WZ $VERSION $CUT
./submit.sh Wj $VERSION $CUT
./submit.sh ZZ $VERSION $CUT

exit
