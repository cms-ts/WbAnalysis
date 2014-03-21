#!/bin/sh

VERSION=v03

if [ ! -z "$1" ]; then
  VERSION=$1
fi

./merge.sh Wj_patgen $VERSION
./merge.sh W1j_patgen $VERSION
./merge.sh W2j_patgen $VERSION
./merge.sh W3j_patgen $VERSION
./merge.sh W4j_patgen $VERSION

./merge.sh W_patgen-all $VERSION

exit
