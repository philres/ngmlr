// To compile:
//   gcc -g -O2 example.c libminimap.a -lz

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "_kseq.h"
#include "minimap.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "Usage: minimap-lite <target.fa> <query.fa>\n");
		return 1;
	}
	
	// open query file for reading; you may use your favorite FASTA/Q parser
	gzFile f = gzopen(argv[2], "r");
	assert(f);
	kseq_t *ks = kseq_init(f);

	// create index for target; we are creating one index for all target sequence
	int n_threads = 4, w = 10, k = 15;
	mm_idx_t *mi = mm_idx_build(argv[1], w, k, n_threads);
	assert(mi);

	// mapping
	mm_mapopt_t opt;
	mm_mapopt_init(&opt); // initialize mapping parameters
	mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
	while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
		const mm_reg1_t *reg;
		int j, n_reg;
		// get all hits for the query
		reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &opt, 0);
		// traverse hits and print them out
		for (j = 0; j < n_reg; ++j) {
			const mm_reg1_t *r = &reg[j];
			printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
			printf("%s\t%d\t%d\t%d\t%d\t%d\n", mi->name[r->rid], mi->len[r->rid], r->rs, r->re, r->len, r->cnt);
		}
	}
	mm_tbuf_destroy(tbuf);

	// deallocate index and close the query file
	mm_idx_destroy(mi);
	kseq_destroy(ks);
	gzclose(f);
	return 0;
}
