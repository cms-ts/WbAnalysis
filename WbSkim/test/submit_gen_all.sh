#!/bin/sh

VERSION=v03

if [ ! -z "$1" ]; then
  VERSION=$1
fi

./submit_gen.sh Wj $VERSION
./submit_gen.sh W1j $VERSION
./submit_gen.sh W2j $VERSION
./submit_gen.sh W3j $VERSION
./submit_gen.sh W4j $VERSION

exit
