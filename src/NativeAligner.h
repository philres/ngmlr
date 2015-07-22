#ifndef __NGMNATIVE_H__
#define __NGMNATIVE_H__

#include "IAlignment.h"
#include "IConfig.h"

#include <string>
#include <emmintrin.h>

typedef int** matrix;

class NativeAligner : public IAlignment
{
protected:
	int m_corr_size;
	int m_corr_size_mem;

	int m_query_len;
	int m_score_match;
	int m_score_mismatch;
	int m_score_gap_open;
	int m_score_gap_extend;

	__m128i* __restrict__ corr_align;
	__m128i* __restrict__ corr_gap_up;

	//Temporary
	matrix scores_align;
	matrix scores_ga;
	matrix scores_gb;
	matrix ptr[3];
	int max_cell_ref_index;
	int max_cell_read_index;

	
public:
	void ScoreAffineLocalSSE(char const * const * const ref, char const * const * const read,float * const results);
	int ScoreAffineLocal(const char* ref, const char* read);
	int ScoreAffineLocalUnopt(const char* ref, const char* read);
	int AlignAffineLocal(const char* ref, const char* read, std::string& aln_ref, std::string& aln_read );

	int ForwardAffineLocal(const char* ref, const char* read );
	void BackwardAffineLocal(const char* ref, const char* read, Align& alignment);

public:
	 NativeAligner();
	~NativeAligner();

	void setCorridorSize(int corr);

	int GetScoreBatchSize() const { return 8192; }
	int GetAlignBatchSize() const { return 8192; }

	int BatchScore(
		int const mode, 
		int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList,
		float * const results,void * extData);

	int BatchAlign(
		int const mode,
		int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList,
		Align * const results,void * extData);

	int SingleAlign(
			int const mode,
			int const corridor,
			char const * const refSeq,
			char const * const qrySeq,
			Align & result,
			void * extData);

	typedef int score_t;

friend class Benchmark;

private:
	void init(int const qryLen, int const corridorLen);
};

#endif
