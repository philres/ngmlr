#!/bin/bash

#set -x

NAME="deterministic mapping quality"
BIN=`ls bin/ | grep -v "debug" | grep "[0-9]$" | sort -n | tail -n 1`
PARAMETER=" --skip-write -x pacbio --no-progress -t 4 -R 0.01 "

echo "Test: $NAME"

for i in {1..5}
do
	bin/${BIN}/ngmlr $PARAMETER -r test/data/test_3/reference.fasta.gz -q test/data/test_3/read.fa.gz 2> /dev/null | samtools view -Sb - 2> /dev/null | bedtools bamtobed | sort > test/data/test_3/found.bed
	diff test/data/test_3/expected.bed test/data/test_3/found.bed || exit 1
	rm test/data/test_3/found.bed
done

echo "Success"
