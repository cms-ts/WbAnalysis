#!/bin/sh

#QUEUE=ts_cms
#QUEUE=short
#QUEUE=normal
QUEUE=normal_io

VERSION=v12
CUT=0

DATADIR=/gpfs/grid/srm/cms/store/user/schizzi/grid

WORKDIR=/gpfs/cms/users/schizzi/Wbb2012/test

OUTDIR=/gpfs/cms/users/schizzi/Wbb2012/test/data

if [ $# -eq 0 ]; then
  echo 'Usage: submit.sh jobdir [version] [cut]'
  exit
fi

JOBDIR=$1

if [ ! -z "$2" ]; then
  VERSION=$2
fi

if [ ! -z "$3" ]; then
  CUT=$3
fi

if [ ! -e $DATADIR/$VERSION/$JOBDIR ]; then
  echo "ERROR: $DATADIR/$VERSION/$JOBDIR does not exist !"
  exit
fi

if [ "$CUT" == "0" ]; then
  rm -fr $OUTDIR/$VERSION/$JOBDIR
  mkdir -p $OUTDIR/$VERSION/$JOBDIR
  cd $OUTDIR/$VERSION/$JOBDIR
else
  rm -fr $OUTDIR/$VERSION.$CUT/$JOBDIR
  mkdir -p $OUTDIR/$VERSION.$CUT/$JOBDIR
  cd $OUTDIR/$VERSION.$CUT/$JOBDIR
fi

find $DATADIR/$VERSION/$JOBDIR -maxdepth 1 -name '*.root' | \
xargs -n 100 \
bsub -q $QUEUE -J $JOBDIR -e /dev/null -o /dev/null $WORKDIR/job.sh $CUT

exit
