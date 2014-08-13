#!/bin/sh

VERSION=v10

if [ ! -z "$1" ]; then
  VERSION=$1
fi

./submit_gen.sh Wj_gen $VERSION
./submit_gen.sh W1j_gen $VERSION
./submit_gen.sh W2j_gen $VERSION
./submit_gen.sh W3j_gen $VERSION
./submit_gen.sh W4j_gen $VERSION

./submit_gen.sh Wbb_gen $VERSION

exit
