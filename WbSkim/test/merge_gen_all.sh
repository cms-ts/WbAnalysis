#!/bin/sh

VERSION=v08

if [ ! -z "$1" ]; then
  VERSION=$1
fi

./merge.sh Wj_gen $VERSION
./merge.sh W1j_gen $VERSION
./merge.sh W2j_gen $VERSION
./merge.sh W3j_gen $VERSION
./merge.sh W4j_gen $VERSION

./merge.sh Wj_gen-all $VERSION

exit
