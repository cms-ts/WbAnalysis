#!/bin/sh

VERSION=v01
USER=$USER

do_check=1
if [ ! -z "$1" ]; then
  if [ "$1" == "-n" ]; then
    do_check=0
    shift
  fi
fi

do_size=0
if [ ! -z "$1" ]; then
  if [ "$1" == "-s" ]; then
    do_size=1
    shift
  fi
fi

if [ ! -z "$1" ]; then
  VERSION=$1
fi

if [ ! -z "$2" ]; then
  USER=$2
fi

DIR=/gpfs/cms/users/schizzi/Wbb2012/test/data/$VERSION/

[ "$USER" == "dellaric" ] && DIR=/gpfs/cms/users/schizzi/Wbb2012/test/GDR/data/$VERSION/

if [ ! -e $DIR ]; then
  echo "ERROR: $DIR does not exist !"
  exit
fi

if [ ! -z "$3" ]; then
  DIRS=$3
else
  DIRS=`ls $DIR | grep -v root`
fi

NT=0

for D in $DIRS; do

  FILES=`find $DIR/$D/ -name job.log`

  N=`find $DIR/$D/ -name job.log | wc -w`

  NT=$((NT+N))

  SIZE=""
  if [ $do_size -eq 1 ]; then
    SIZE=`du -sh $DIR/$D/ | awk '{print $1}'`
  fi

  printf "%30s\t%s%6i%10s\n" $D ":" $N $SIZE

  for F in $FILES; do

    if [ -e $F ]; then

      if [ ! -s $F ]; then
        echo "ERROR: empty file "
        ls -lt $F
      fi

      if [ $do_check -eq 1 ]; then
        E=`grep "cmsRun: command not found" $F`
        if [ ! -z "$E" ]; then
          echo "ERROR: cmsRun: command not found in "$F
        fi
        E=`grep "Begin Fatal Exception" $F`
        if [ ! -z "$E" ]; then
          echo "ERROR: exception in "$F
        fi
        E=`grep "fatal system signal" $F`
        if [ ! -z "$E" ]; then
          echo "ERROR: fatal system signal in "$F
        fi
        E=`grep "Abort" $F`
        if [ ! -z "$E" ]; then
          echo "ERROR: abort in "$F
        fi
      fi

    fi

  done

done

SIZE=""
if [ $do_size -eq 1 ]; then
  if [ ! -z "$3" ]; then
    SIZE=`du -sh $DIR/$3 | awk '{print $1}'`
  else
    SIZE=`du -sh $DIR | awk '{print $1}'`
  fi
fi

printf "%30s\t%s%6i%10s\n" "TOTAL" ":" $NT $SIZE

exit
