/*
 * EndToEndAffine.h
 *
 *  Created on: Oct 22, 2013
 *      Author: philipp_
 */

#ifndef ENDTOENDAFFINE_H_
#define ENDTOENDAFFINE_H_

#include "IAlignment.h"
#include "ILog.h"
#include "IConfig.h"

#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/file.h>  // for printint strings

using seqan::Dna5String;
using seqan::Score;
using seqan::Simple;
using seqan::AlignConfig;
using seqan::StringSet;
using seqan::length;

class EndToEndAffine: public IAlignment {

    typedef seqan::String<char> TSequence;
    typedef seqan::StringSet<TSequence> TStringSet;
    typedef seqan::Gaps<TSequence, seqan::ArrayGaps> TGaps;
    typedef seqan::Iterator<TGaps>::Type TGapsIterator;
    typedef seqan::Iterator<seqan::String<int> >::Type TIterator;

public:
	EndToEndAffine() : m_AlignMode(Config.GetInt(MODE, 0, 1)), lDiag(0), uDiag(Config.GetInt("corridor")) {

		scoringScheme = Score<float, Simple>(Config.GetFloat(MATCH_BONUS), Config.GetFloat(MISMATCH_PENALTY) * -1.0f, Config.GetFloat(GAP_EXTEND_PENALTY) * -1.0f, Config.GetFloat(GAP_READ_PENALTY) * -1.0f);

		//scoringScheme = Score<float, Simple>(5.0, -4.0, -0.5, -10.0);
		//scoringScheme = Score<float, Simple>(10, -15, -5, -33);
	}
	virtual ~EndToEndAffine() {}

	virtual int GetScoreBatchSize() const {
		return 1024;
	}
	virtual int GetAlignBatchSize() const {
		return 1024;
	}

	virtual int BatchScore(int const mode, int const batchSize, char const * const * const refSeqList, char const * const * const qrySeqList, float * const results, void * extData);

	virtual int BatchAlign(int const mode, int const batchSize, char const * const * const refSeqList, char const * const * const qrySeqList, Align * const results, void * extData);

	void convertToCIGAR(TGaps gapsText, TGaps gapsPattern, Align & result, int const readLen);

private:
	bool const m_AlignMode;
	int const lDiag;
	int const uDiag;
	Score<float, Simple> scoringScheme;
	seqan::AlignConfig<true, true, false, true> alignConfig;

};

#endif /* ENDTOENDAFFINE_H_ */
