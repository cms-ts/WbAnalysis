#!/bin/sh

VERSION=v10
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

./submit.sh TTbar_FullLept $VERSION $CUT
./submit.sh TTbar_SemiLept $VERSION $CUT

./submit.sh T_s $VERSION $CUT
./submit.sh T_t $VERSION $CUT
./submit.sh T_tW $VERSION $CUT
./submit.sh TBar_s $VERSION $CUT
./submit.sh TBar_t $VERSION $CUT
./submit.sh TBar_tW $VERSION $CUT

./submit.sh Wj $VERSION $CUT
./submit.sh W1j $VERSION $CUT
./submit.sh W2j $VERSION $CUT
./submit.sh W3j $VERSION $CUT
./submit.sh W4j $VERSION $CUT

./submit.sh Wbb $VERSION $CUT

./submit.sh WW $VERSION $CUT
./submit.sh WZ $VERSION $CUT
./submit.sh ZZ $VERSION $CUT

exit
