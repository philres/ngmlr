#!/bin/bash

NAME="simple read length test"
BIN="ngmlr-"`grep -o "set( NGM_VERSION_MAJOR [0-9]* )" CMakeLists.txt | cut -d " " -f 3`"."`grep -o "set( NGM_VERSION_MINOR [0-9]* )" CMakeLists.txt | cut -d " " -f 3`"."`grep -o "set( NGM_VERSION_BUILD [^ ]* )" CMakeLists.txt | grep -v debug | cut -d " " -f 3`
PARAMETER=" --skip-write "

echo "Test: $NAME ($BIN)"

bin/${BIN}/ngmlr $PARAMETER -r test/data/test_2/ref_chr21_20kb.fa -q test/data/test_2/reads_100_2200bp.fa 2> /dev/null | samtools view -Sb - 2> /dev/null | bedtools bamtobed > test/data/test_2/found.bed

diff test/data/test_2/expected.bed test/data/test_2/found.bed  && echo "Success"

rm test/data/test_2/found.bed