#!/bin/bash

#set -x

NAME="read name length"
BIN="ngmlr-"`grep -o "set( NGM_VERSION_MAJOR [0-9]* )" CMakeLists.txt | cut -d " " -f 3`"."`grep -o "set( NGM_VERSION_MINOR [0-9]* )" CMakeLists.txt | cut -d " " -f 3`"."`grep -o "set( NGM_VERSION_BUILD [^ ]* )" CMakeLists.txt | grep -v debug | cut -d " " -f 3`
PARAMETER=" --skip-write "

echo "Test: $NAME ($BIN)"

bin/${BIN}/ngmlr $PARAMETER -r test/data/test_1/ref_chr6_140kb.fa -q test/data/test_1/long_name.fa 2> /dev/null | grep -v "^@" | cut -f 1 > test/data/test_1/found.txt

diff test/data/test_1/expected.txt test/data/test_1/found.txt && echo "Success"

rm test/data/test_1/found.txt