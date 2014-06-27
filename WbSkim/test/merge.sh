#!/bin/sh

VERSION=v09

if [ $# -eq 0 ]; then
  echo 'Usage: merge.sh jobdir [version]'
  exit
fi

JOBDIR=$1

if [ ! -z "$2" ]; then
  VERSION=$2
fi

WORKDIR=/gpfs/cms/users/schizzi/Wbb2012/test/data

if [ ! -e $WORKDIR/$VERSION ]; then
  echo 'ERROR: version "'$VERSION'" does not exist !'
  exit
fi

export ROOT_HIST=0

unset PYTHIA8175DATA
unset G4NEUTRONXS

cd $CMS_PATH/slc6_amd64_gcc472/cms/cmssw/CMSSW_6_2_7
eval `scramv1 runtime -sh`
cd - > /dev/null

opts="-T -v 0"

if [ -d $WORKDIR/$VERSION/$JOBDIR ]; then
  echo 'Preparing '$WORKDIR/$VERSION/$JOBDIR.root
  rm -f $WORKDIR/$VERSION/$JOBDIR.root
  hadd $opts $WORKDIR/$VERSION/$JOBDIR.root $WORKDIR/$VERSION/$JOBDIR/LSFJOB_*/rootTuple_*.root
elif [ "$JOBDIR" == "data-all" ]; then
  echo 'Preparing '$WORKDIR/$VERSION/SingleElectron_2012_merge.root
  rm -f $WORKDIR/$VERSION/SingleElectron_2012_merge.root
  hadd $opts $WORKDIR/$VERSION/SingleElectron_2012_merge.root $WORKDIR/$VERSION/SingleElectron_2012*.root
  echo 'Preparing '$WORKDIR/$VERSION/SingleMu_2012_merge.root
  rm -f $WORKDIR/$VERSION/SingleMu_2012_merge.root
  hadd $opts $WORKDIR/$VERSION/SingleMu_2012_merge.root $WORKDIR/$VERSION/SingleMu*_2012*.root
elif [ "$JOBDIR" == "T-all" ]; then
  rm -f $WORKDIR/$VERSION/T_merge.root
  root -l -q -b hmerge.C\(\"$WORKDIR\",\"$VERSION\",\"T\"\)
elif [ "$JOBDIR" == "TTbar-all" ]; then
  rm -f $WORKDIR/$VERSION/TTbar_merge.root
  root -l -q -b hmerge.C\(\"$WORKDIR\",\"$VERSION\",\"TTbar\"\)
elif [ "$JOBDIR" == "W-all" ]; then
  rm -f $WORKDIR/$VERSION/W_merge.root
  root -l -q -b hmerge.C\(\"$WORKDIR\",\"$VERSION\",\"W\"\)
elif [ "$JOBDIR" == "Wj_patgen-all" ]; then
  rm -f $WORKDIR/$VERSION/W_patgen_merge.root
  root -l -q -b hmerge.C\(\"$WORKDIR\",\"$VERSION\",\"W_patgen\"\)
elif [ "$JOBDIR" == "Wj_gen-all" ]; then
  rm -f $WORKDIR/$VERSION/W_gen_merge.root
  root -l -q -b hmerge.C\(\"$WORKDIR\",\"$VERSION\",\"W_gen\"\)
else
  echo 'ERROR: jobdir "'$JOBDIR'" does not exist !'
fi

exit
