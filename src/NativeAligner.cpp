#include "NativeAligner.h"
#include "NGM.h"

#include <iostream>
#include <stdio.h>
#include <string>
using namespace std;

#define max(a,b) (a>b?a:b)
#define min(a,b) (a<b?a:b)

int STOP = -1;
int ALIGN = 0;
int GA = 1;
int GB = 2;

matrix make_matrix(int size, int init_val = 0) {
	matrix mat = new int*[size];
	for (int i = 0; i < size; ++i) {
		mat[i] = new int[size];
		for (int j = 0; j < size; ++j)
			mat[i][j] = init_val;
	}

	return mat;
}

void delete_matrix(matrix mat, int size) {
	for (int i = 0; i < size; ++i) {
		delete[] mat[i];
	}
	delete[] mat;
}

void debug_matrix(int len, int** scores, int** ptr, std::string ref,
		std::string read, std::string title = "", int xmark = -1,
		int ymark = -1) {

	std::cout << title << std::endl;

	std::string ref_n = "-" + ref;
	std::string read_n = "-" + read;

	std::cout << "  ";
	for (int i = 0; i < len; i++) {
		std::cout << "      " << ref_n[i];
	}
	std::cout << " (< read) (v ref)";

	std::cout << endl;

	for (int j = 0; j < len; j++) {
		std::cout << read_n[j] << "  ";

		for (int i = 0; i < len; ++i) {
			char ptr_char;
			if (ptr[i][j] == STOP)
				ptr_char = '!';
			else if (ptr[i][j] == ALIGN)
				ptr_char = '\\';
			else if (ptr[i][j] == GA)
				ptr_char = '<';
			else if (ptr[i][j] == GB)
				ptr_char = '^';

			if (i == xmark && j == ymark)
				ptr_char = '*';

			printf(" %4d%c ", scores[i][j], ptr_char);
		}
		std::cout << endl;
		std::cout << endl;
	}
}

#include <stdlib.h>
#define malloc_align(ptr,size,align)  posix_memalign(((void * *)&ptr), align, size)

void NativeAligner::init(int const qryLen, int const corridorLen) {
	m_corr_size = 0;
	m_query_len = qryLen + 5;
	setCorridorSize(corridorLen);
	//Temporary, for simple version of backtracking/CIGAR computation
	scores_align = make_matrix(m_query_len + 1);
	scores_ga = make_matrix(m_query_len + 1, 0);
	scores_gb = make_matrix(m_query_len + 1, 0);
	ptr[ALIGN] = make_matrix(m_query_len + 1, STOP);
	ptr[GA] = make_matrix(m_query_len + 1);
	ptr[GB] = make_matrix(m_query_len + 1);
}

NativeAligner::NativeAligner() {
	Log.Message("Using NativeAligner (SSE)");

	m_score_match = Config.GetInt("match_bonus");
	m_score_mismatch = -Config.GetInt("mismatch_penalty");
	m_score_gap_open = -Config.GetInt("gap_read_penalty");
	m_score_gap_extend = -Config.GetInt("gap_extend_penalty"); //Fixme

	corr_align = 0;
	corr_gap_up = 0;

//	init(Config.GetInt("qry_max_len"), Config.GetInt("corridor"));
}


void NativeAligner::setCorridorSize(int corr) {
	m_corr_size = corr;
	m_corr_size_mem = corr + 1;

	free(corr_align);
	free(corr_gap_up);

	//corr_align = new score_t[m_corr_size+1];
	malloc_align(corr_align, m_corr_size_mem * sizeof(__m128i ), 16);
	for (int i = 0; i < m_corr_size_mem; i++)
		corr_align[i] = _mm_setzero_si128();

	malloc_align(corr_gap_up, m_corr_size_mem * sizeof(__m128i ), 16);
	for (int i = 0; i < m_corr_size_mem; i++)
		corr_gap_up[i] = _mm_setzero_si128();
}

NativeAligner::~NativeAligner() {

}

int NativeAligner::BatchScore(int const mode, int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList, float * const results,
		void * extData) {
	int i = 0;
	for (; i < batchSize - 7; i += 8) {
		ScoreAffineLocalSSE(refSeqList + i, qrySeqList + i, results + i);
	}

	if (i > batchSize - 7) {
		i = batchSize - 8;
		ScoreAffineLocalSSE(refSeqList + i, qrySeqList + i, results + i);
	}

	/*for( int i = 0; i < batchSize; i ++ )
	 {
	 results[i] = ScoreAffineLocal( *(refSeqList + i),*(qrySeqList + i));
	 }*/

	return batchSize;

}
#include<stdlib.h>
int NativeAligner::BatchAlign(int const mode, int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList, Align * const results,
		void * extData) {

	for (int i = 0; i < batchSize; ++i) {
		Align& aln = results[i];
		aln.Score = ForwardAffineLocal(refSeqList[i], qrySeqList[i]);
		BackwardAffineLocal(refSeqList[i], qrySeqList[i], aln);
	}

	return batchSize;
}

int NativeAligner::SingleAlign(int const mode, int const corridor,
		char const * const refSeq, char const * const qrySeq, Align & result,
		void * extData) {

	int qryLen = strlen(qrySeq);

	init(qryLen, corridor);

	Align& aln = result;
	aln.Score = ForwardAffineLocal(refSeq, qrySeq);
	BackwardAffineLocal(refSeq, qrySeq, aln);

	return 1;
}

int NativeAligner::AlignAffineLocal(const char* ref, const char* read,
		std::string& aln_ref, std::string& aln_read) {
	int len = m_query_len;

	int max_cell_score = 0;
	int max_cell_read_index = 0;
	int max_cell_ref_index = 0;

	for (int read_index = 1; read_index < len + 1; read_index++) {
		for (int ref_index = 1; ref_index < len + 1; ref_index++) {
			int score_gap_left_open = scores_align[read_index - 1][ref_index]
					+ m_score_gap_open;
			int score_gap_left_extend = scores_ga[read_index - 1][ref_index]
					+ m_score_gap_extend;

			if (score_gap_left_open >= score_gap_left_extend)
				ptr[GA][read_index][ref_index] = ALIGN;
			else
				ptr[GA][read_index][ref_index] = GA;

			scores_ga[read_index][ref_index] = max(score_gap_left_open,
					score_gap_left_extend);

			///////////////////////////////////////////////////////

			int score_gap_up_open = scores_align[read_index][ref_index - 1]
					+ m_score_gap_open;
			int score_gap_up_extend = scores_gb[read_index][ref_index - 1]
					+ m_score_gap_extend;

			if (score_gap_up_open >= score_gap_up_extend)
				ptr[GB][read_index][ref_index] = ALIGN;
			else
				ptr[GB][read_index][ref_index] = GB;

			scores_gb[read_index][ref_index] = max(score_gap_up_open,
					score_gap_up_extend);

			///////////////////////////////////////////////////////

			int score_align = scores_align[read_index - 1][ref_index - 1]
					+ ((read[read_index - 1] == ref[ref_index - 1]) ?
							m_score_match : m_score_mismatch);

			//std::cout << "read " << read[read_index-1] << " vs ref "  << ref[ref_index-1] << ": " << scores_align[read_index-1][ref_index-1] << " score add " << ((read[read_index-1] == ref[ref_index-1]) ? m_score_match : m_score_mismatch) << endl;
			//std::cout << score_align << endl;

			int score_cell = max(score_align, 0);
			score_cell = max(scores_ga[read_index][ref_index], score_cell);
			score_cell = max(scores_gb[read_index][ref_index], score_cell);
			scores_align[read_index][ref_index] = score_cell;

			if (score_cell == score_align)
				ptr[ALIGN][read_index][ref_index] = ALIGN;
			else if (score_cell == scores_ga[read_index][ref_index])
				ptr[ALIGN][read_index][ref_index] = GA;
			else if (score_cell == scores_gb[read_index][ref_index])
				ptr[ALIGN][read_index][ref_index] = GB;

			if (score_cell == 0)
				ptr[ALIGN][read_index][ref_index] = STOP;

			//std::cout << "ptr " << read_index << "," << ref_index << " setto " << ptr_align[read_index][ref_index] << endl;

			if (score_cell > max_cell_score) {
				max_cell_score = score_cell;
				max_cell_read_index = read_index;
				max_cell_ref_index = ref_index;
				//std::cout << "(" << max_cell_read_index << "," << max_cell_ref_index << ")" << endl;
				//std::cout << "(" << read[read_index-1] << "," << ref[ref_index-1] << ")" << endl;
			}

		}

		//std::cout<<endl;

	}

	//std::cout << "ptr " <<1 << "," << 1 << " eq " << ptr_align[1][1] << endl;

	// std::cout << "ALN" << std::endl;
	// debug_matrix(len+1,scores_align,ptr_align,ref,read);

	// std::cout << "GA=LEFT" << std::endl;
	// debug_matrix(len+1,scores_ga,ptr_ga,ref,read);

	// std::cout << "GB=UP" << std::endl;
	// debug_matrix(len+1,scores_gb,ptr_gb,ref,read);

	int score = max_cell_score;
	int curr_read_index = max_cell_read_index;
	int curr_ref_index = max_cell_ref_index;

	int** ptr_matrix = ptr[ALIGN];

	int aln_read_index = max_cell_read_index - 1;
	int aln_ref_index = max_cell_ref_index - 1;

	while (true) {
		int action = ptr_matrix[curr_read_index][curr_ref_index];
		//std::cout << "(" << curr_read_index << "," << curr_ref_index << ") act=" << action << " ";
		//cout << scores_align[curr_read_index][curr_ref_index];

		if (ptr_matrix == ptr[ALIGN]) {
			if (action == GA || action == GB) {
				ptr_matrix = ptr[action];
				continue;
			}

			if (action == STOP)
				break;

			aln_ref.push_back(ref[aln_ref_index--]);
			aln_read.push_back(read[aln_read_index--]);
			curr_read_index--;
			curr_ref_index--;
			//debug_matrix(len+1,scores_align,ptr[ALIGN],ref,read,"ALN",curr_read_index,curr_ref_index);
		} else if (ptr_matrix == ptr[GA]) {
			aln_read.push_back(read[aln_read_index--]);
			aln_ref.push_back('-');
			curr_read_index--;
			//debug_matrix(len+1,scores_ga,ptr[GA],ref,read,"GA/LEFT",curr_read_index,curr_ref_index);
		} else if (ptr_matrix == ptr[GB]) {
			aln_ref.push_back(ref[aln_ref_index--]);
			aln_read.push_back('-');
			curr_ref_index--;
			//debug_matrix(len+1,scores_gb,ptr[GB],ref,read,"GB/UP",curr_read_index,curr_ref_index);
		}

		ptr_matrix = ptr[action];

		//cout << endl;
		//getchar();
	}

	//debug_matrix(len+1,scores_align,ptr[ALIGN],ref,read,"ALN",curr_read_index,curr_ref_index);

	//std::cout << "max_cell: " << max_cell_score << std::endl;

	//delete_matrix(scores_align,len+1);
	//delete_matrix(scores_ga,len+1);
	//delete_matrix(scores_gb,len+1);
	//delete_matrix(ptr[ALIGN],len+1);
	//delete_matrix(ptr[GA],len+1);
	//delete_matrix(ptr[GB],len+1);

	return max_cell_score;

}

int NativeAligner::ScoreAffineLocalUnopt(const char* ref, const char* read) {
	int corr_array_size = m_corr_size;
	int* corr_align = new score_t[corr_array_size];
	for (int i = 0; i < corr_array_size; i++)
		corr_align[i] = 0;
	int* corr_gb = new score_t[corr_array_size];
	for (int i = 0; i < corr_array_size; i++)
		corr_gb[i] = 0;

	int max_score = 0;
	int ref_offset = -corr_array_size / 2;

	for (int read_index = 0; read[read_index]; read_index++) {
		int score_left = 0;
		int score_left_gap = 0;
		int score_up_gap = 0;

		for (int ref_index = 0; ref_index < corr_array_size; ref_index++) {
			if (ref_offset + ref_index < 0)
				continue;
			if (ref[ref_offset + ref_index] == 0)
				break;

			score_left_gap = max(score_left_gap + m_score_gap_extend, 0);
			score_left_gap = max(score_left + m_score_gap_open, score_left_gap);

			if (ref_index + 1 < corr_array_size) {
				score_up_gap = max(corr_gb[ref_index + 1] + m_score_gap_extend,
						0);
				score_up_gap = max(corr_align[ref_index + 1] + m_score_gap_open,
						score_up_gap);
			} else {
				score_up_gap = 0;
			}

			score_left = max(score_left_gap, 0);
			score_left = max(score_up_gap, score_left);

			int score_align =
					(ref[ref_offset + ref_index] == read[read_index]) ?
							m_score_match : m_score_mismatch;
			score_left = max(corr_align[ref_index] + score_align, score_left);

			corr_align[ref_index] = score_left;
			corr_gb[ref_index] = score_up_gap;

			//std::cout << "Pos(" << read_index << "," << ref_offset + ref_index << ") : (" << read[read_index] << "," << ref[ref_offset+ref_index] << ")= " << score_left << " (alnscore " << score_align << ")" << std::endl;

			max_score = max(max_score, score_left);
		}

		ref_offset++;
	}

	delete[] corr_align;
	delete[] corr_gb;
	return max_score;
}

int NativeAligner::ScoreAffineLocal(const char* ref, const char* read) {
	int len = m_query_len;

	int corr_array_size = m_corr_size;
	int* corr_align = new score_t[corr_array_size + 1];
	for (int i = 0; i < corr_array_size + 1; i++)
		corr_align[i] = 0;
	int* corr_gap_up = new score_t[corr_array_size + 1];
	for (int i = 0; i < corr_array_size + 1; i++)
		corr_gap_up[i] = 0;

	int max_score = 0;
	int ref_offset = -corr_array_size >> 1;

	for (int read_index = 0; read_index < len; read_index++) {
		int score_left = 0;
		int score_left_gap = 0;
		int score_up_gap = 0;

		int ref_index_start = max(0, -ref_offset);
		int ref_index_max = min(corr_array_size, len - ref_offset);

		for (int ref_index = ref_index_start; ref_index < ref_index_max;
				ref_index++) {
			score_left_gap = max(score_left_gap + m_score_gap_extend, 0);
			score_left_gap = max(score_left + m_score_gap_open, score_left_gap);

			score_up_gap = max(corr_gap_up[ref_index + 1] + m_score_gap_extend,
					0);
			score_up_gap = max(corr_align[ref_index + 1] + m_score_gap_open,
					score_up_gap);

			score_left = max(score_left_gap, 0);
			score_left = max(score_up_gap, score_left);

			int score_align =
					(ref[ref_offset + ref_index] == read[read_index]) ?
							m_score_match : m_score_mismatch;
			score_left = max(corr_align[ref_index] + score_align, score_left);

			corr_align[ref_index] = score_left;
			corr_gap_up[ref_index] = score_up_gap;

			max_score = max(max_score, score_left);
		}

		ref_offset++;
	}

	delete[] corr_align;
	delete[] corr_gap_up;
	return max_score;
}

void NativeAligner::ScoreAffineLocalSSE(char const * const * const ref_,
		char const * const * const read_, float* results) {
	for (int i = 0; i < m_corr_size_mem; i++)
		corr_align[i] = _mm_setzero_si128();
	for (int i = 0; i < m_corr_size_mem; i++)
		corr_gap_up[i] = _mm_setzero_si128();

	int len = m_query_len;

	const __m128i score_gap_extend = _mm_set1_epi16(m_score_gap_extend);
	const __m128i score_gap_open = _mm_set1_epi16(m_score_gap_open);
	const __m128i score_match = _mm_set1_epi16(m_score_match);
	const __m128i score_mismatch = _mm_set1_epi16(m_score_mismatch);
	const __m128i zero = _mm_setzero_si128();

	__m128i max_score = _mm_setzero_si128();
	int ref_offset = -m_corr_size >> 1;

	for (int read_index = 0; read_index < len; read_index++) {
		const __m128i read_curr = _mm_set_epi16(read_[0][read_index],
				read_[1][read_index], read_[2][read_index],
				read_[3][read_index], read_[4][read_index],
				read_[5][read_index], read_[6][read_index],
				read_[7][read_index]);

		__m128i score_left = _mm_setzero_si128();
		__m128i score_left_gap = _mm_setzero_si128();
		__m128i score_up_gap = _mm_setzero_si128();

		int ref_index_start = max(0, -ref_offset);
		int ref_index_max = min(m_corr_size, len - ref_offset);

		for (int ref_index = ref_index_start; ref_index < ref_index_max;
				ref_index++) {
			int ref_actual_index = ref_offset + ref_index;
			const __m128i ref_curr = _mm_set_epi16(ref_[0][ref_actual_index],
					ref_[1][ref_actual_index], ref_[2][ref_actual_index],
					ref_[3][ref_actual_index], ref_[4][ref_actual_index],
					ref_[5][ref_actual_index], ref_[6][ref_actual_index],
					ref_[7][ref_actual_index]);

			score_left_gap = _mm_max_epi16(
					_mm_add_epi16(score_left_gap, score_gap_extend), zero);
			score_left_gap = _mm_max_epi16(
					_mm_add_epi16(score_left, score_gap_open), score_left_gap);

			score_up_gap = _mm_max_epi16(
					_mm_add_epi16(corr_gap_up[ref_index + 1], score_gap_extend),
					zero);
			score_up_gap = _mm_max_epi16(
					_mm_add_epi16(corr_align[ref_index + 1], score_gap_open),
					score_up_gap);

			score_left = _mm_max_epi16(score_left_gap, zero);
			score_left = _mm_max_epi16(score_up_gap, score_left);

			__m128i cmp = _mm_cmpeq_epi16(read_curr, ref_curr);
			__m128i score_align = _mm_or_si128(_mm_and_si128(cmp, score_match),
					_mm_andnot_si128(cmp, score_mismatch));

			score_left = _mm_max_epi16(
					_mm_add_epi16(corr_align[ref_index], score_align),
					score_left);

			corr_align[ref_index] = score_left;
			corr_gap_up[ref_index] = score_up_gap;

			max_score = _mm_max_epi16(max_score, score_left);
		}

		ref_offset++;
	}

	short* scores_convert = (short*) &max_score;
	results[0] = scores_convert[7];
	results[1] = scores_convert[6];
	results[2] = scores_convert[5];
	results[3] = scores_convert[4];
	results[4] = scores_convert[3];
	results[5] = scores_convert[2];
	results[6] = scores_convert[1];
	results[7] = scores_convert[0];
}

int NativeAligner::ForwardAffineLocal(const char* ref, const char* read) {
	for (int i = 0; i < m_query_len; ++i) {
		for (int j = 0; j < m_query_len; ++j) {
			scores_align[i][j] = 0;
			ptr[GA][i][j] = 0;
			ptr[GB][i][j] = 0;
		}
	}

	int len = m_query_len;
	int max_cell_score = 0;

	//int ref_offset = -m_corr_size >> 1;

	for (int read_index = 1; read_index < len + 1; read_index++) {
		//int ref_index_start = max(1,-ref_offset);
		//int ref_index_max = min(m_corr_size,(len+1)-ref_offset);

		//for( int ref_index = ref_index_start; ref_index < ref_index_max; ref_index ++ )
		for (int ref_index = 1; ref_index < len + 1; ref_index++) {
			int actual_ref_index = ref_index; // + ref_offset;

			int score_gap_left_open = scores_align[read_index - 1][ref_index]
					+ m_score_gap_open;
			int score_gap_left_extend = scores_ga[read_index - 1][ref_index]
					+ m_score_gap_extend;

			if (score_gap_left_open >= score_gap_left_extend)
				ptr[GA][read_index][ref_index] = ALIGN;
			else
				ptr[GA][read_index][ref_index] = GA;

			scores_ga[read_index][ref_index] = max(score_gap_left_open,
					score_gap_left_extend);

			///////////////////////////////////////////////////////

			int score_gap_up_open = scores_align[read_index][ref_index - 1]
					+ m_score_gap_open;
			int score_gap_up_extend = scores_gb[read_index][ref_index - 1]
					+ m_score_gap_extend;

			if (score_gap_up_open >= score_gap_up_extend)
				ptr[GB][read_index][ref_index] = ALIGN;
			else
				ptr[GB][read_index][ref_index] = GB;

			scores_gb[read_index][ref_index] = max(score_gap_up_open,
					score_gap_up_extend);

			///////////////////////////////////////////////////////

			int score_align = scores_align[read_index - 1][ref_index - 1]
					+ ((read[read_index - 1] == ref[actual_ref_index - 1]) ?
							m_score_match : m_score_mismatch);

			int score_cell = max(score_align, 0);
			score_cell = max(scores_ga[read_index][ref_index], score_cell);
			score_cell = max(scores_gb[read_index][ref_index], score_cell);
			scores_align[read_index][ref_index] = score_cell;

			if (score_cell == score_align)
				ptr[ALIGN][read_index][ref_index] = ALIGN;
			else if (score_cell == scores_ga[read_index][ref_index])
				ptr[ALIGN][read_index][ref_index] = GA;
			else if (score_cell == scores_gb[read_index][ref_index])
				ptr[ALIGN][read_index][ref_index] = GB;

			if (score_cell == 0)
				ptr[ALIGN][read_index][ref_index] = STOP;

			if (score_cell > max_cell_score) {
				max_cell_score = score_cell;
				max_cell_read_index = read_index;
				max_cell_ref_index = ref_index;
			}
		}

		//ref_offset ++;
	}
	return max_cell_score;
}

void NativeAligner::BackwardAffineLocal(const char* ref, const char* read,
		Align& aln) {
	std::string cigar_raw;
	std::string md_raw;

	int curr_read_index = max_cell_read_index;
	int curr_ref_index = max_cell_ref_index;

	int** ptr_matrix = ptr[ALIGN];

	int aln_read_index = max_cell_read_index - 1;
	int aln_ref_index = max_cell_ref_index - 1;
	int aln_index = max_cell_read_index - 1;

	bool in_deletion = false;

	while (true) {
		//if(dbg) std::cin.get();

		int action = ptr_matrix[curr_read_index][curr_ref_index];

		if (ptr_matrix == ptr[ALIGN]) {
			if (action == GA || action == GB) {
				ptr_matrix = ptr[action];
				continue;
			}

			if (action == STOP) {
				break;
			}

			char base_ref = ref[aln_ref_index--];
			char base_read = read[aln_read_index--];
			if (base_ref == base_read) {
				aln.Identity++;
				md_raw = "M" + md_raw;
			} else {
				aln.NM++;

				std::string helper;
				helper.push_back(base_ref);
				md_raw = helper + md_raw;
			}

			cigar_raw = "M" + cigar_raw;
			curr_read_index--;
			curr_ref_index--;
			aln_index--;
			in_deletion = false;

		} else if (ptr_matrix == ptr[GA]) {
			cigar_raw = "I" + cigar_raw;
			curr_read_index--;
			aln_read_index--;
			aln.NM++;

			in_deletion = false;

		} else if (ptr_matrix == ptr[GB]) {
			cigar_raw = "D" + cigar_raw;

			if (!in_deletion)
				md_raw = "^" + md_raw;
			in_deletion = true;

			curr_ref_index--;
			aln_ref_index--;
			aln_index--;
			aln.NM++;
		}

		ptr_matrix = ptr[action];
	}

	aln.QStart = 0;
	aln.QEnd = 0;

	std::string cigar;
	std::string md;
	int clip_front = curr_read_index;
	for (int i = 0; i < clip_front; ++i) {
		cigar_raw = "S" + cigar_raw;
		//aln.QStart ++;
	}

	for (int i = max_cell_read_index; read[i] != 0 && i < m_query_len; ++i) {
		cigar_raw = cigar_raw + "S";
		//aln.QEnd ++;	
	}

	int raw_len = cigar_raw.size();
	for (int i = 0; i < cigar_raw.size(); ++i)
		if (cigar_raw[i] == 'D')
			raw_len--;

	char current = cigar_raw[0];
	int current_start = 0;
	for (int i = 1; i < cigar_raw.size(); ++i) {
		if (current != cigar_raw[i]) {
			int count = i - current_start;

			char count_str[255];
			sprintf((char*) count_str, "%d", count);
			cigar += std::string(count_str);
			cigar.push_back(current);

			current = cigar_raw[i];
			current_start = i;
		}
	}

	int count = cigar_raw.size() - current_start;

	if (count > 0) {
		char count_str[255];
		sprintf((char*) count_str, "%d", count);
		cigar += std::string(count_str);
		cigar.push_back(current);
	}

	//MD finalization
	{
		char current = md_raw[0];
		int current_start = 0;
		for (int i = 1; i < md_raw.size(); ++i) {
			if (md_raw[i] == 'M') {
				if (current != 'M')
					current_start = i - 1;
				current = 'M';
			} else {
				if (current == 'M') {
					int count = i - current_start;

					char count_str[255];
					sprintf((char*) count_str, "%d", count);
					md += std::string(count_str);
				}

				md.push_back(md_raw[i]);
				current = md_raw[i];
			}
		}

		int count = md_raw.size() - current_start;

		if (current == 'M' && count > 0) {
			char count_str[255];
			sprintf((char*) count_str, "%d", count);
			md += std::string(count_str);
		}
	}

	strcpy(aln.pBuffer1, cigar.c_str());
	strcpy(aln.pBuffer2, md.c_str());

	/*std::cout<<"Read:";
	 for( int i =0; i < strlen(read); i ++ )
	 {
	 std::cout<<read[i];
	 }
	 std::cout << std::endl << "Ref: ";
	 for( int i =0; i < strlen(read); i ++ )
	 {
	 std::cout<<ref[i];
	 }
	 std::cout << std::endl;
	 std::cout << "Score: " << aln.Score << std::endl;
	 std::cout << "Readlen: " << strlen(read) << std::endl;
	 std::cout << "CIGAR: " << cigar << std::endl;
	 debug_matrix(m_query_len+1,scores_align,ptr[ALIGN],read,ref,"",max_cell_read_index,max_cell_ref_index);*/
}
