/*
 * cout-reads.h
 *
 *  Created on: Apr 4, 2014
 *      Author: fritz
 */

#ifndef COUT_READS_H_
#define COUT_READS_H_

int count_reads(int argc, char **argv);
int NextSAMRead(IParser * parser, int const id,SAMRecord * & read);

#endif /* COUT_READS_H_ */
