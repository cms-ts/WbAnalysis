#!/bin/sh

#QUEUE=ts_cms
#QUEUE=short
#QUEUE=normal
QUEUE=normal_io

VERSION=v08

DATADIR=/gpfs/grid/srm/cms/store/user/schizzi/grid

WORKDIR=/gpfs/cms/users/schizzi/Wbb2012/test

OUTDIR=/gpfs/cms/users/schizzi/Wbb2012/test/data

if [ $# -eq 0 ]; then
  echo 'Usage: submit.sh jobdir [version]'
  exit
fi

JOBDIR=$1

if [ ! -z "$2" ]; then
  VERSION=$2
fi

if [ ! -e $DATADIR/$VERSION/$JOBDIR ]; then
  echo "ERROR: $DATADIR/$VERSION/$JOBDIR does not exist !"
  exit
fi

rm -fr $OUTDIR/$VERSION/${JOBDIR}_patgen
mkdir -p $OUTDIR/$VERSION/${JOBDIR}_patgen
cd $OUTDIR/$VERSION/${JOBDIR}_patgen

find $DATADIR/$VERSION/$JOBDIR -maxdepth 1 -name '*.root' | \
xargs -n 50 \
bsub -q $QUEUE -J $JOBDIR -e /dev/null -o /dev/null $WORKDIR/job_patgen.sh

exit
