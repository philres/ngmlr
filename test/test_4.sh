#!/bin/bash

#set -x

NAME="primary alignment test"
BIN=`ls bin/ | grep -v "debug" | grep "[0-9]$" | sort -n | tail -n 1`
PARAMETER=" --skip-write -x pacbio --no-progress -t 4 "

echo "Test: $NAME"

bin/${BIN}/ngmlr $PARAMETER -r test/data/test_4/reference.fasta.gz -q test/data/test_4/read.fa.gz 2> /dev/null | samtools view -S - 2> /dev/null | sort | cut -f 1,2,3,4,5 > test/data/test_4/found.txt

diff test/data/test_4/expected.txt test/data/test_4/found.txt || exit 1

rm test/data/test_4/found.txt

echo "Success"
