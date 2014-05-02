/*
 * Options.h
 *
 *  Created on: Jul 7, 2012
 *      Author: philipp_
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <getopt.h>

static char const * getopt_short = "c:o:q:r:t:gs:m:f:k:pleI:X:i:n:R:C:b1:2:d:Q:";

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
		{ MODE, 						required_argument, 0, 'm' },
		{ "kmer", 						required_argument, 0, 'k' },
		{ "min-insert-size",			required_argument, 0, 'I' },
		{ "max-insert-size",			required_argument, 0, 'X' },
		{ "paired", 					no_argument,       0, 'p' },
		{ "topn", 						required_argument, 0, 'n' },
		{ "kmer-skip", 					required_argument, 0, 0 },
		{ "affine", 					no_argument, 0, 0 },
		{ "max-consec-indels",		required_argument, 0, 'C' },
		{ "pre-only", 					no_argument,       0, 0 },
		{ "skip-env", 					no_argument,       0, 0 },
		{ "bam",    					no_argument,       0, 'b' },
		{ "color",    					no_argument,       0, 0 },
		{ "max-equal", 					required_argument, 0, 0 },
		{ "search-table-length", 		required_argument, 0, 0 },
		{ "local",				 		no_argument, 0, 'l' },
		{ "end-to-end",			 		no_argument, 0, 'e' },
		{ "bs-mapping", 				no_argument		 , 0, 0 },
		{ "min-identity",               required_argument, 0, 'i'},
		{ "min-residues",               required_argument, 0, 'R'},
		{ "min-score",					required_argument, 0, 'S'},
		{ "min-mq",  					required_argument, 0, 'Q'},
		{ "parse-all",                  no_argument      , 0, 0 },
		{ "hard-clip",                  no_argument      , 0, 0 },
		{ "silent-clip",                no_argument      , 0, 0 },
		{ "bs-cutoff",                  required_argument, 0, 0 },
		{ "kmer-min",                   required_argument, 0, 0 },
		{ "match-bonus",				required_argument, 0, 0 },
		{ "match-bonus-tt",				required_argument, 0, 0 },
		{ "match-bonus-tc",				required_argument, 0, 0 },
		{ "mismatch-penalty",			required_argument, 0, 0 },
		{ "gap-read-penalty",			required_argument, 0, 0 },
		{ "gap-ref-penalty",			required_argument, 0, 0 },
		{ "gap-extend-penalty",			required_argument, 0, 0 },
		{ "max-cmrs",      		        required_argument, 0, 0 },
		{ "fast-pairing",  		        no_argument,       0, 0 },
		{ "no-unal",	  		        no_argument,       0, 0 },
		{ "skip-save",  		        no_argument,       0, 0 },
		{ "no-progress",  		        no_argument,       0, 0 },
		{ "pair-score-cutoff",	        required_argument, 0, 0 },
		{ "strata",      		        no_argument,       0, 0 },
		{ "pe-delimiter",  		        required_argument, 0, 'd' },
		{ "keep-tags",      	        no_argument,       0, 0 },
		{ "skip-mate-check",	        no_argument,       0, 0 },
		{ "rg-id",	       				required_argument, 0, 0 },
		{ "rg-cn",       				required_argument, 0, 0 },
		{ "rg-ds",      				required_argument, 0, 0 },
		{ "rg-dt",      				required_argument, 0, 0 },
		{ "rg-fo",     	    			required_argument, 0, 0 },
		{ "rg-ks",     	 	    		required_argument, 0, 0 },
		{ "rg-lb",     			    	required_argument, 0, 0 },
		{ "rg-pg",     				    required_argument, 0, 0 },
		{ "rg-pi",      				required_argument, 0, 0 },
		{ "rg-pl",      				required_argument, 0, 0 },
		{ "rg-pu",      				required_argument, 0, 0 },
		{ "rg-sm",      				required_argument, 0, 0 },
#ifdef DEBUGLOG
		{ "log",        				required_argument, 0, 0 },
		{ "log-lvl",      				required_argument, 0, 0 },
#endif
		{ "argos",   	  		        no_argument,       0, 0 },
		{ "argos-min-score",			required_argument, 0, 0 },
	0 };

#endif /* OPTIONS_H_ */
