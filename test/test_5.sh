#!/bin/bash

#set -x

NAME="max query name length test"
BIN="ngmlr-"`grep -o "set( NGM_VERSION_MAJOR [0-9]* )" CMakeLists.txt | cut -d " " -f 3`"."`grep -o "set( NGM_VERSION_MINOR [0-9]* )" CMakeLists.txt | cut -d " " -f 3`"."`grep -o "set( NGM_VERSION_BUILD [0-9]* )" CMakeLists.txt | cut -d " " -f 3`
PARAMETER=" --skip-write -x pacbio --no-progress -t 4 "

echo "Test: $NAME ($BIN)"

LEN=`bin/${BIN}/ngmlr $PARAMETER -r test/data/test_5/reference.fasta.gz -q test/data/test_5/read.fa.gz 2> /dev/null | grep -v "^@" | cut -f 1 | wc -c`
echo $LEN

if [ "$LEN" -gt 254 ]
then
	exit 1		
fi

echo "Success"
