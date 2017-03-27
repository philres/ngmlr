#!/bin/bash

NAME="simple read length test"
BIN=`ls bin/ | grep -v "debug" | grep "[0-9]$" | sort -n | tail -n 1`
PARAMETER=" --skip-write "

echo "Test: $NAME"

bin/${BIN}/ngmlr $PARAMETER -r test/data/test_2/ref_chr21_20kb.fa -q test/data/test_2/reads_100_2200bp.fa 2> /dev/null | samtools view -Sb - 2> /dev/null | bedtools bamtobed > test/data/test_2/found.bed

diff test/data/test_2/expected.bed test/data/test_2/found.bed  && echo "Success"

