/**
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Contact: philipp.rescheneder@univie.ac.at
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
	EndToEndAffine() : lDiag(0), uDiag(Config.getReadPartCorridor()) {

		scoringScheme = Score<float, Simple>(Config.getScoreMatch(), Config.getScoreMismatch(), Config.getScoreExtend(), Config.getScoreGapOpen());

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

	virtual int SingleScore(int const mode, int const corridor, char const * const refSeq, char const * const qrySeq, float & result, void * extData);

	virtual int SingleAlign(int const mode, int const corridor, char const * const refSeq, char const * const qrySeq, Align & result, void * extData);

	void convertToCIGAR(int const mode, TGaps gapsText, TGaps gapsPattern, Align & result, int const readLen);

private:
	int const lDiag;
	int const uDiag;
	Score<float, Simple> scoringScheme;
	seqan::AlignConfig<true, true, false, true> alignConfig;

};

#endif /* ENDTOENDAFFINE_H_ */
