/*
 * Options.h
 *
 *  Created on: Jul 7, 2012
 *      Author: philipp_
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <getopt.h>

char const * getopt_short = "c:o:q:r:d:t:gs:m:f:k:pI:X:i:n:R:C:b1:2:";

//Here '-' has to be used and not '_'. When querying the Config the '-' has to be replaced by '_'.
//When passing parameters trough the config file '-' and '_' can both be used.
static const struct option long_options[] =
	{
		{ "config", 					required_argument, 0, 'c' },
		{ "output", 					required_argument, 0, 'o' },
		{ "qry",    					required_argument, 0, 'q' },
		{ "qry1",	    				required_argument, 0, '1' },
		{ "qry2", 	   					required_argument, 0, '2' },
		{ "ref", 						required_argument, 0, 'r' },
		{ "cpu-threads", 				required_argument, 0, 't' },
		{ "gpu", 						no_argument      , 0, 'g' },
		{ "sensitivity", 				required_argument, 0, 's' },
		{ "mode", 						required_argument, 0, 'm' },
		{ "kmer", 						required_argument, 0, 'k' },
		{ "min-insert-size",			required_argument, 0, 'I' },
		{ "max-insert-size",			required_argument, 0, 'X' },
		{ "paired", 					no_argument,       0, 'p' },
		{ "topn", 						required_argument, 0, 'n' },
		{ "kmer-skip", 					required_argument, 0, 0 },
		{ "max-consecutive-indels",		required_argument, 0, 'C' },
		{ "pre-only", 					no_argument,       0, 0 },
		{ "skip-env", 					no_argument,       0, 0 },
		{ "bam",    					no_argument,       0, 'b' },
		{ "color",    					no_argument,       0, 0 },
		{ "max-equal", 					required_argument, 0, 0 },
		{ "search-table-length", 		required_argument, 0, 0 },
		{ "bs-mapping", 				no_argument		 , 0, 0 },
		{ "min-identity",               required_argument, 0, 'i'},
		{ "min-residues",               required_argument, 0, 'R'},
		{ "parse-all",                  no_argument      , 0, 0 },
		{ "hard-clip",                  no_argument      , 0, 0 },
		{ "silent-clip",                no_argument      , 0, 0 },
		{ "bs-cutoff",                  required_argument, 0, 0 },
		{ "kmer-min",                   required_argument, 0, 0 },
		{ "score-match",				required_argument, 0, 0 },
		{ "score-match-tt",				required_argument, 0, 0 },
		{ "score-match-tc",				required_argument, 0, 0 },
		{ "score-mismatch",				required_argument, 0, 0 },
		{ "score-gap-read",				required_argument, 0, 0 },
		{ "score-gap-ref",				required_argument, 0, 0 },
		{ "max-cmrs",      		        required_argument, 0, 0 },
		{ "fast-pairing",  		        no_argument,       0, 0 },
		{ "no-unal",	  		        no_argument,       0, 0 },
		{ "skip-save",  		        no_argument,       0, 0 },
		{ "no-progress",  		        no_argument,       0, 0 },
		{ "pair-score-cutoff",	        required_argument, 0, 0 },
		{ "strata",      		        no_argument,       0, 0 },
		{ "pe-delimiter",  		        required_argument, 0, 'd' },
	0 };

#endif /* OPTIONS_H_ */
