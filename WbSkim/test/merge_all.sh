#!/bin/sh

VERSION=v09

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

./merge.sh TTbar_FullLept $VERSION $CUT
./merge.sh TTbar_SemiLept $VERSION $CUT

./merge.sh TTbar-all $VERSION

./merge.sh T_s
./merge.sh T_t
./merge.sh T_tW
./merge.sh TBar_s
./merge.sh TBar_t
./merge.sh TBar_tW

./merge.sh T-all $VERSION

./merge.sh Wj $VERSION
./merge.sh W1j $VERSION
./merge.sh W2j $VERSION
./merge.sh W3j $VERSION
./merge.sh W4j $VERSION

./merge.sh W-all $VERSION

./merge.sh Wbb $VERSION

./merge.sh WW $VERSION
./merge.sh WZ $VERSION
./merge.sh ZZ $VERSION

exit
