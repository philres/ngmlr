#!/bin/bash

#set -x

NAME="reads from stdin"
BIN="ngmlr-"`grep -o "set( NGM_VERSION_MAJOR [0-9]* )" CMakeLists.txt | cut -d " " -f 3`"."`grep -o "set( NGM_VERSION_MINOR [0-9]* )" CMakeLists.txt | cut -d " " -f 3`"."`grep -o "set( NGM_VERSION_BUILD [^ ]* )" CMakeLists.txt | grep -v debug | cut -d " " -f 3`
PARAMETER=" --skip-write -x pacbio --no-progress -t 4 "

echo "Test: $NAME ($BIN)"

LEN=`cat test/data/test_6/read.fa.gz | gunzip | bin/${BIN}/ngmlr $PARAMETER -r test/data/test_6/reference.fasta.gz 2> /dev/null | samtools view -Sc - 2> /dev/null`

if [ "$LEN" -eq "5" ]
then
	echo "Success"
else
	exit 1		
fi

