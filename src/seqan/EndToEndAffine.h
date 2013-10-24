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
	EndToEndAffine() {}
	virtual ~EndToEndAffine() {}

	virtual int GetScoreBatchSize() const {
		return 1024;
	}
	virtual int GetAlignBatchSize() const {
		return 1024;
	}

	virtual int BatchScore(int const mode, int const batchSize, char const * const * const refSeqList, char const * const * const qrySeqList, float * const results, void * extData) {
		Score<int, Simple> scoringScheme(10, -15, -5, -33);
		seqan::AlignConfig<true, true, false, true> alignConfig;
		int lDiag = 0;
		static int uDiag = Config.GetInt("corridor");

		for(int i = 0; i < batchSize; ++i) {
			//Dna5String seqH = "CGATT";
			//Dna5String seqV = "CGAAATT";
			results[i] = (float)globalAlignmentScore(TSequence(refSeqList[i]), TSequence(qrySeqList[i]), scoringScheme, alignConfig, lDiag, uDiag);
			//results[i] = (float)globalAlignmentScore(TSequence(refSeqList[i]), TSequence(qrySeqList[i]), seqan::MyersBitVector());
		}
		return batchSize;
	}

	void convertToCIGAR(TGaps gapsText, TGaps gapsPattern, Align & result) {

		TGapsIterator itGapsPattern = seqan::begin(gapsPattern);
		TGapsIterator itGapsEnd = seqan::end(gapsPattern);

		int match = 0;
		int mismatch = 0;
		int total = 0;

		// Remove trailing gaps in pattern.
		int count = 0;
		while(isGap(--itGapsEnd))
			++count;
		seqan::setClippedEndPosition(gapsPattern, length(gapsPattern) - count);

		// Remove leading gaps in pattern.
		if(isGap(itGapsPattern))
		{
			seqan::setClippedBeginPosition(gapsPattern, countGaps(itGapsPattern));
			int offest = countGaps(itGapsPattern);
			seqan::setClippedBeginPosition(gapsText, offest);
			result.PositionOffset = offest;
		}

		// Reinitilaize the iterators.
		TGapsIterator itGapsText = seqan::begin(gapsText);
		itGapsPattern = seqan::begin(gapsPattern);
		itGapsEnd = seqan::end(gapsPattern);

		// Use a stringstream to construct the cigar string.
		std::stringstream cigar;
		int numChar = 0;
		while (itGapsPattern != itGapsEnd)
		{
			// Count insertions.
			if (isGap(itGapsText))
			{
				int numGaps = countGaps(itGapsText);
				cigar << numGaps << "I";
				itGapsText += numGaps;
				itGapsPattern += numGaps;
				total += numGaps;
				continue;
			}
			// Count deletions.
			if (isGap(itGapsPattern))
			{
				int numGaps = countGaps(itGapsPattern);
				cigar << numGaps << "D";
				itGapsText += numGaps;
				itGapsPattern += numGaps;
				total += numGaps;
				continue;
			}
			// Count matches.
			//while (*itGapsText == *itGapsPattern && itGapsPattern != itGapsEnd)
			while (!isGap(itGapsText) && !isGap(itGapsPattern) && itGapsPattern != itGapsEnd)
			{
				if(*itGapsText == *itGapsPattern) {
					match += 1;
				} else {
					mismatch += 1;
				}
				total += 1;
				++numChar;
				++itGapsText;
				++itGapsPattern;
			}
			if (numChar != 0)
			{
				cigar << numChar << "M";
				numChar = 0;
				continue;
			}
			//// Count mismatches.
			//while (*itGapsText != *itGapsPattern && itGapsPattern != itGapsEnd)
			//{
			//	++numChar;
			//	++itGapsText;
			//	++itGapsPattern;
			//}
			//if (numChar != 0)
			//	cigar << numChar << "S";
			//numChar = 0;
		}
		// Output the hit position in the text, the total number of edits and the corresponding cigar string.
		//::std::cout << "Hit at position  " << *it << "\ttotal edits: " << abs(score) << "\tcigar: " << cigar.str() << ::std::endl;
		//std::cout << "--------------------------------------------------------------------------------------------------" << std::endl;

		result.QStart = 0;
		result.QEnd = 0;
		result.Identity = match * 1.0f / total;
		result.NM = mismatch;
		strcpy(result.pBuffer1, cigar.str().c_str());


	}

	virtual int BatchAlign(int const mode, int const batchSize, char const * const * const refSeqList, char const * const * const qrySeqList, Align * const results, void * extData) {
		int lDiag = 0;
		static int uDiag = Config.GetInt("corridor");

		for(int i = 0; i < batchSize; ++i) {

			TGaps gapsText;
			TGaps gapsPattern;

			//Log.Message("Ref:  %s", refSeqList[i]);
			//Log.Message("Read: %s", qrySeqList[i]);
			assignSource(gapsText, TSequence(refSeqList[i]));
			assignSource(gapsPattern, TSequence(qrySeqList[i]));

			Score<int, Simple> scoringScheme(10, -15, -5, -33);
			seqan::AlignConfig<true, true, false, true> alignConfig;

			float result = (float)globalAlignment(gapsText, gapsPattern, scoringScheme, alignConfig, lDiag, uDiag);
			//int result = (float)localAlignment(gapsText, gapsPattern, scoringScheme, lDiag, uDiag);

			//std::cout << "1: " << gapsText << std::endl << "2:" << gapsPattern << std::endl;
			//Log.Message("Score: %f", result);

			//std::cerr << "Score: " << result << std::endl;
			//std::cerr << "The resulting alignment is\n"
			//<< convertToCIGAR(gapsText, gapsPattern) << std::endl;

			convertToCIGAR(gapsText, gapsPattern, results[i]);

		}
		return batchSize;
	}

};

#endif /* ENDTOENDAFFINE_H_ */
