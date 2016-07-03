#ifndef __BITPALALIGNER_H__
#define __BITPALALIGNER_H__

#include "IAlignment.h"
#include "IConfig.h"

class BitpalAligner: public IAlignment {

protected:

	int Score(char const * ref, char const * read);

public:
	BitpalAligner();
	~BitpalAligner();

//	void setCorridorSize(int corr);

	int GetScoreBatchSize() const {
		return 8192;
	}
	int GetAlignBatchSize() const {
		return 8192;
	}

	int BatchScore(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList, float * const results,
			void * extData);

	int BatchAlign(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList, Align * const results,
			void * extData);

	int SingleAlign(int const mode, int const corridor,
			char const * const refSeq, char const * const qrySeq,
			Align & result, void * extData);

};

#endif//__BITPALALIGNER_H__
