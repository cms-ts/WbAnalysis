#!/bin/sh

VERSION=v01

if [ ! -z "$1" ]; then
  VERSION=$1
fi

./merge.sh SingleElectron_2012A_22Jan13 $VERSION
./merge.sh SingleElectron_2012B_22Jan13 $VERSION
./merge.sh SingleElectron_2012C_22Jan13 $VERSION
./merge.sh SingleElectron_2012D_22Jan13 $VERSION

./merge.sh SingleMu_2012A_22Jan13 $VERSION
./merge.sh SingleMu_2012B_22Jan13 $VERSION
./merge.sh SingleMu_2012C_22Jan13 $VERSION
./merge.sh SingleMu_2012D_22Jan13 $VERSION

./merge.sh data-all $VERSION

./merge.sh DYJetsToLL $VERSION
./merge.sh TTbar $VERSION
./merge.sh WZ $VERSION
./merge.sh Wj $VERSION
./merge.sh W1j $VERSION
./merge.sh W2j $VERSION
./merge.sh W3j $VERSION
./merge.sh W4j $VERSION
./merge.sh WW $VERSION
./merge.sh ZZ $VERSION

exit
