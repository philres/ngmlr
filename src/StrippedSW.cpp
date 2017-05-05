/**
 * Contact: philipp.rescheneder@gmail.com
 */

#include "StrippedSW.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "ssw.h"

//	Print the BLAST like output.
static void ssw_write(const s_align* a, const char* ref_seq,
		const char* read_seq, const int8_t* table) {

	fprintf(stdout,
			"optimal_alignment_score: %d\tsub-optimal_alignment_score: %d\t",
			a->score1, a->score2);
	if (a->ref_begin1 + 1)
		fprintf(stdout, "target_begin: %d\t", a->ref_begin1 + 1);
	fprintf(stdout, "target_end: %d\t", a->ref_end1 + 1);
	if (a->read_begin1 + 1)
		fprintf(stdout, "query_begin: %d\t", a->read_begin1 + 1);
	fprintf(stdout, "query_end: %d\n\n", a->read_end1 + 1);
	if (a->cigar) {
		int32_t c = 0, left = 0, e = 0, qb = a->ref_begin1, pb = a->read_begin1;
		uint32_t i;
		while (e < a->cigarLen || left > 0) {
			int32_t count = 0;
			int32_t q = qb;
			int32_t p = pb;
			fprintf(stdout, "Target: %8d    ", q + 1);
			for (c = e; c < a->cigarLen; ++c) {
				char letter = cigar_int_to_op(a->cigar[c]);
				uint32_t length = cigar_int_to_len(a->cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left : length;
				for (i = 0; i < l; ++i) {
					if (letter == 'I')
						fprintf(stdout, "-");
					else {
						fprintf(stdout, "%c", *(ref_seq + q));
						++q;
					}
					++count;
					if (count == 60)
						goto step2;
				}
			}
			step2: fprintf(stdout, "    %d\n                    ", q);
			q = qb;
			count = 0;
			for (c = e; c < a->cigarLen; ++c) {
				char letter = cigar_int_to_op(a->cigar[c]);
				uint32_t length = cigar_int_to_len(a->cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left : length;
				for (i = 0; i < l; ++i) {
					if (letter == 'M') {
						if (table[(int) *(ref_seq + q)]
								== table[(int) *(read_seq + p)])
							fprintf(stdout, "|");
						else
							fprintf(stdout, "*");
						++q;
						++p;
					} else {
						fprintf(stdout, "*");
						if (letter == 'I')
							++p;
						else
							++q;
					}
					++count;
					if (count == 60) {
						qb = q;
						goto step3;
					}
				}
			}
			step3: p = pb;
			fprintf(stdout, "\nQuery:  %8d    ", p + 1);
			count = 0;
			for (c = e; c < a->cigarLen; ++c) {
				char letter = cigar_int_to_op(a->cigar[c]);
				uint32_t length = cigar_int_to_len(a->cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left : length;
				for (i = 0; i < l; ++i) {
					if (letter == 'D')
						fprintf(stdout, "-");
					else {
						fprintf(stdout, "%c", *(read_seq + p));
						++p;
					}
					++count;
					if (count == 60) {
						pb = p;
						left = l - i - 1;
						e = (left == 0) ? (c + 1) : c;
						goto end;
					}
				}
			}
			e = c;
			left = 0;
			end: fprintf(stdout, "    %d\n\n", p);
		}
	}
}

/* This table is used to transform nucleotide letters into numbers. */
static const int8_t nt_table[128] = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 0,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

int StrippedSW::BatchScore(int const mode, int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList, float * const results,
		void * extData) {
	for (int i = 0; i < batchSize; ++i) {

		char const * const ref_seq = refSeqList[i];
		char const * const read_seq = qrySeqList[i];

		//Log.Message("Ref:  %s", ref_seq);
		//Log.Message("Read: %s", read_seq);

		int read_len = strlen(read_seq) + 1;
		int ref_len = strlen(ref_seq) + 1;

		if (read_len >= maxSeqLen || ref_len >= maxSeqLen) {
			results[i] = -1.0f;
		} else {
			s_profile* profile;
			s_align* result;

			for (int32_t m = 0; m < read_len; ++m)
				num[m] = nt_table[(int) read_seq[m]];
			num[read_len] = nt_table[(int) '\0'];

			profile = ssw_init(num, read_len, mat, 5, 1);
			for (int32_t m = 0; m < ref_len; ++m)
				ref_num[m] = nt_table[(int) ref_seq[m]];
			ref_num[ref_len] = nt_table[(int) '\0'];

			// Only the 8 bit of the flag is setted. ssw_align will always return the best alignment beginning position and cigar.
			result = ssw_align(profile, ref_num, ref_len, gap_open, gap_extension, 0, 0, 0, 0);
			//ssw_write(result, ref_seq, read_seq, nt_table);
			//fprintf(stderr, "%d\n", result->score1);
			results[i] = result->score1;

			align_destroy(result);
			init_destroy(profile);
		}

	}
	return batchSize;
}

int StrippedSW::SingleScore(int const mode, int const corridor,
		char const * const refSeq, char const * const qrySeq, float & resultScore,
		void * extData) {
	char const * const ref_seq = refSeq;
	char const * const read_seq = qrySeq;

	//Log.Message("Ref:  %s", ref_seq);
	//Log.Message("Read: %s", read_seq);

	int read_len = strlen(read_seq) + 1;
	int ref_len = strlen(ref_seq) + 1;

	if(read_len >= maxSeqLen || ref_len >= maxSeqLen) {
		resultScore = -1.0f;
		return 0;
	}

	s_profile* profile;
	s_align * result;

	for (int32_t m = 0; m < read_len; ++m)
		num[m] = nt_table[(int) read_seq[m]];
	num[read_len] = nt_table[(int) '\0'];

	profile = ssw_init(num, read_len, mat, 5, 1);
	for (int32_t m = 0; m < ref_len; ++m)
		ref_num[m] = nt_table[(int) ref_seq[m]];
	ref_num[ref_len] = nt_table[(int) '\0'];

	// Only the 8 bit of the flag is setted. ssw_align will always return the best alignment beginning position and cigar.
	result = ssw_align(profile, ref_num, ref_len, gap_open, gap_extension, 0, 0,
			0, 0);
	//ssw_write(result, ref_seq, read_seq, nt_table);
	//fprintf(stderr, "%d\n", result->score1);
	resultScore = result->score1;

	align_destroy(result);
	init_destroy(profile);

	return 1;
}

int StrippedSW::SingleAlign(int const mode, int const corridor,
		char const * const refSeq, char const * const qrySeq, Align & results,
		void * extData) {
	char const * const ref_seq = refSeq;
	char const * const read_seq = qrySeq;

	Log.Message("Aligning with StrippedSW");
//	Log.Message("Ref:  %s", ref_seq);
//	Log.Message("Read: %s", read_seq);

	int read_len = strlen(read_seq) + 1;
	int ref_len = strlen(ref_seq) + 1;

	s_profile* profile;
	s_align * result;

	for (int32_t m = 0; m < read_len; ++m)
		num[m] = nt_table[(int) read_seq[m]];
	num[read_len] = nt_table[(int) '\0'];

	profile = ssw_init(num, read_len, mat, 5, 1);
	for (int32_t m = 0; m < ref_len; ++m)
		ref_num[m] = nt_table[(int) ref_seq[m]];
	ref_num[ref_len] = nt_table[(int) '\0'];

	// Only the 8 bit of the flag is setted. ssw_align will always return the best alignment beginning position and cigar.
	result = ssw_align(profile, ref_num, ref_len, gap_open, gap_extension, 1, 0,
			0, 0);
	//ssw_write(result, ref_seq, read_seq, nt_table);
	//Log.Message("%d", result->score1);
	Align & align = results;
	char * outcigar = align.pBuffer1;

	read_len -= 1;
	int sum = 0;
	align.QStart = result->read_begin1;
	sum += align.QStart;
	if (align.QStart > 0) {
		outcigar += sprintf(outcigar, "%d%c", align.QStart, 'S');
	}
	for (int i = 0; i < result->cigarLen; ++i) {
		if (cigar_int_to_op(result->cigar[i]) != 'D') {
			sum += cigar_int_to_len(result->cigar[i]);
		}
		outcigar += sprintf(outcigar, "%d%c",
				cigar_int_to_len(result->cigar[i]),
				cigar_int_to_op(result->cigar[i]));
	}
	align.QEnd = (read_len - result->read_end1 - 1);
	sum += align.QEnd;
	if (align.QEnd > 0) {
		outcigar += sprintf(outcigar, "%d%c", align.QEnd, 'S');
	}
	*outcigar = '\0';

	if (sum != read_len) {
		fprintf(stderr, "%d == %d -- %d %d %d\n", sum, read_len, align.QStart,
				align.QEnd, read_len);
		ssw_write(result, ref_seq, read_seq, nt_table);
		fprintf(stderr, "%s\n", align.pBuffer1);
	}

	//results[i] = result->score1;
	align.PositionOffset = result->ref_begin1;
	align.Identity = 1.0f;
	align.NM = 0;

	align_destroy(result);
	init_destroy(profile);

	return 1;
}

int StrippedSW::BatchAlign(int const mode, int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList, Align * const results,
		void * extData) {
	for (int i = 0; i < batchSize; ++i) {

		char const * const ref_seq = refSeqList[i];
		char const * const read_seq = qrySeqList[i];

		//Log.Message("Ref:  %s", ref_seq);
		//Log.Message("Read: %s", read_seq);

		int read_len = strlen(read_seq) + 1;
		int ref_len = strlen(ref_seq) + 1;

		s_profile* profile;
		s_align* result;

		for (int32_t m = 0; m < read_len; ++m)
			num[m] = nt_table[(int) read_seq[m]];
		num[read_len] = nt_table[(int) '\0'];

		profile = ssw_init(num, read_len, mat, 5, 1);
		for (int32_t m = 0; m < ref_len; ++m)
			ref_num[m] = nt_table[(int) ref_seq[m]];
		ref_num[ref_len] = nt_table[(int) '\0'];

		// Only the 8 bit of the flag is setted. ssw_align will always return the best alignment beginning position and cigar.
		result = ssw_align(profile, ref_num, ref_len, gap_open, gap_extension,
				1, 0, 0, 0);
		//ssw_write(result, ref_seq, read_seq, nt_table);
		//Log.Message("%d", result->score1);
		Align & align = results[i];
		char * outcigar = align.pBuffer1;

		read_len -= 1;
		int sum = 0;
		align.QStart = result->read_begin1;
		sum += align.QStart;
		if (align.QStart > 0) {
			outcigar += sprintf(outcigar, "%d%c", align.QStart, 'S');
		}
		for (int i = 0; i < result->cigarLen; ++i) {
			if (cigar_int_to_op(result->cigar[i]) != 'D') {
				sum += cigar_int_to_len(result->cigar[i]);
			}
			outcigar += sprintf(outcigar, "%d%c",
					cigar_int_to_len(result->cigar[i]),
					cigar_int_to_op(result->cigar[i]));
		}
		align.QEnd = (read_len - result->read_end1 - 1);
		sum += align.QEnd;
		if (align.QEnd > 0) {
			outcigar += sprintf(outcigar, "%d%c", align.QEnd, 'S');
		}
		*outcigar = '\0';

		if (sum != read_len) {
			fprintf(stderr, "%d == %d -- %d %d %d\n", sum, read_len,
					align.QStart, align.QEnd, read_len);
			ssw_write(result, ref_seq, read_seq, nt_table);
			fprintf(stderr, "%s\n", align.pBuffer1);
		}

		//results[i] = result->score1;
		align.PositionOffset = result->ref_begin1;
		align.Identity = 1.0f;
		align.NM = 0;

		align_destroy(result);
		init_destroy(profile);

	}
	return batchSize;
}
