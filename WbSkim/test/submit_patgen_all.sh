#!/bin/sh

VERSION=v03

if [ ! -z "$1" ]; then
  VERSION=$1
fi

./submit_patgen.sh Wj $VERSION
./submit_patgen.sh W1j $VERSION
./submit_patgen.sh W2j $VERSION
./submit_patgen.sh W3j $VERSION
./submit_patgen.sh W4j $VERSION

exit
