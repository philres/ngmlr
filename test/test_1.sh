#!/bin/bash

#set -x

NAME="read name length"
BIN=`ls bin/ | grep -v "debug" | grep "[0-9]$" | sort -n | tail -n 1`
PARAMETER=" --skip-write "

echo "Test: $NAME"

bin/${BIN}/ngmlr $PARAMETER -r test/data/test_1/ref_chr6_140kb.fa -q test/data/test_1/long_name.fa 2> /dev/null | grep -v "^@" | cut -f 1 > test/data/test_1/found.txt

diff test/data/test_1/expected.txt test/data/test_1/found.txt && echo "Success"

