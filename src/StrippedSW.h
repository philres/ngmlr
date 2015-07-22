/*
 * EndToEndAffine.h
 *
 *  Created on: Oct 22, 2013
 *      Author: philipp_
 */

#ifndef STRIPPEDSW_H_
#define STRIPPEDSW

#include "ssw.h"

#include "IAlignment.h"
#include "ILog.h"
#include "IConfig.h"


class StrippedSW: public IAlignment {


public:
  StrippedSW() : m_AlignMode(Config.GetInt(MODE, 0, 1)), lDiag(0), uDiag(Config.GetInt("corridor")) {

		//scoringScheme = Score<float, Simple>(Config.GetFloat(MATCH_BONUS), Config.GetFloat(MISMATCH_PENALTY) * -1.0f, Config.GetFloat(GAP_EXTEND_PENALTY) * -1.0f, Config.GetFloat(GAP_READ_PENALTY) * -1.0f);
    gap_open = 3;
    gap_extension = 1;
    match = 2;
    mismatch = 2;
    int32_t l, m, k;	// default parameters for genome sequence alignment
    // initialize scoring matrix for genome sequences
    //  A  C  G  T	N (or other ambiguous code)
    //  2 -2 -2 -2 	0	A
    // -2  2 -2 -2 	0	C
    // -2 -2  2 -2 	0	G
    // -2 -2 -2  2 	0	T
    //	0  0  0  0  0	N (or other ambiguous code)
    mat = (int8_t*)calloc(25, sizeof(int8_t));
    for (l = k = 0; l < 4; ++l) {
      for (m = 0; m < 4; ++m) mat[k++] = l == m ? match : - mismatch;	/* weight_match : -weight_mismatch */
      mat[k++] = 0; // ambiguous base: no penalty
    }
    for (m = 0; m < 5; ++m) mat[k++] = 0;

    num = (int8_t*)malloc(200);	// the read sequence represented in numbers
    ref_num = (int8_t*)malloc(200);	// the read sequence represented in numbers

	}
	virtual ~StrippedSW() {
    free(ref_num);
    free(num);
    free(mat);
  }

	virtual int GetScoreBatchSize() const {
		return 1024;
	}
	virtual int GetAlignBatchSize() const {
		return 1024;
	}

	virtual int BatchScore(int const mode, int const batchSize, char const * const * const refSeqList, char const * const * const qrySeqList, float * const results, void * extData);

	virtual int BatchAlign(int const mode, int const batchSize, char const * const * const refSeqList, char const * const * const qrySeqList, Align * const results, void * extData);

private:
	bool const m_AlignMode;
	int const lDiag;
	int const uDiag;

  int8_t* mat;

  int32_t match, mismatch, gap_open, gap_extension;

  int8_t* num; 
  int8_t* ref_num;
};

#endif /* STRIPPEDSW_H_ */
