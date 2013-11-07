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
	EndToEndAffine() : m_AlignMode(Config.GetInt("mode", 0, 1)), lDiag(0), uDiag(Config.GetInt("corridor")) {

		scoringScheme = Score<float, Simple>(Config.GetFloat("score_match"), Config.GetFloat("score_mismatch"), Config.GetFloat("score_gap_extend"), Config.GetFloat("score_gap_read"));

		//scoringScheme = Score<float, Simple>(5.0, -4.0, -0.5, -10.0);
		//scoringScheme = Score<float, Simple>(10, -15, -5, -33);


		//		Score<int, Simple> scoringScheme(8, -20, -1, -33);
		//Score<int, Simple> scoringScheme(30, -35, -1, -60);
	}
	virtual ~EndToEndAffine() {}

	virtual int GetScoreBatchSize() const {
		return 1024;
	}
	virtual int GetAlignBatchSize() const {
		return 1024;
	}

	virtual int BatchScore(int const mode, int const batchSize, char const * const * const refSeqList, char const * const * const qrySeqList, float * const results, void * extData) {
		for(int i = 0; i < batchSize; ++i) {
			//Dna5String seqH = "CGATT";
			//Dna5String seqV = "CGAAATT";
			switch(m_AlignMode) {
				case 0:
				{
					TGaps gapsText;
					TGaps gapsPattern;
					assignSource(gapsText, TSequence(refSeqList[i]));
					assignSource(gapsPattern, TSequence(qrySeqList[i]));

					results[i] = (float)localAlignment(gapsText, gapsPattern, scoringScheme, lDiag, uDiag);
					//results[i] = (float)globalAlignmentScore(TSequence(refSeqList[i]), TSequence(qrySeqList[i]), scoringScheme, alignConfig, lDiag, uDiag);
				}break;
				case 1:
					results[i] = (float)globalAlignmentScore(TSequence(refSeqList[i]), TSequence(qrySeqList[i]), scoringScheme, alignConfig, lDiag, uDiag);
					//result = (float)globalAlignment(gapsText, gapsPattern, scoringScheme, alignConfig, lDiag, uDiag);
					break;
			}
			//results[i] = (float)globalAlignmentScore(TSequence(refSeqList[i]), TSequence(qrySeqList[i]), seqan::MyersBitVector());
		}
		return batchSize;
	}

	void convertToCIGAR(TGaps gapsText, TGaps gapsPattern, Align & result, int const readLen) {

		result.QStart = clippedBeginPosition(gapsPattern);
		result.QEnd = 0;
		result.PositionOffset = 0;

		TGapsIterator itGapsPattern = seqan::begin(gapsPattern);
		TGapsIterator itGapsEnd = seqan::end(gapsPattern);

		int match = 0;
		int mismatch = 0;
		int total = 0;

		int count = 0;
		if(m_AlignMode == 1) {
			// Remove trailing gaps in pattern.
			while(isGap(--itGapsEnd))
				++count;
			seqan::setClippedEndPosition(gapsPattern, length(gapsPattern) - count);
		}

		// Remove leading gaps in pattern.
		if(isGap(itGapsPattern))
		{
			seqan::setClippedBeginPosition(gapsPattern, countGaps(itGapsPattern));
			int offest = countGaps(itGapsPattern);
			seqan::setClippedBeginPosition(gapsText, offest);
			result.PositionOffset = offest;
		} else {
			result.PositionOffset += clippedBeginPosition(gapsText);
		}

		// Reinitilaize the iterators.
		TGapsIterator itGapsText = seqan::begin(gapsText);
		itGapsPattern = seqan::begin(gapsPattern);
		itGapsEnd = seqan::end(gapsPattern);

		// Use a stringstream to construct the cigar string.
		std::stringstream cigar;
		int numChar = 0;
		if(result.QStart > 0) {
			cigar << result.QStart << "S";
		}
		int numDels = 0;
		while (itGapsPattern != itGapsEnd)
		{
			// Count insertions.
			if (isGap(itGapsText))
			{
				int numGaps = countGaps(itGapsText);
				cigar << numGaps << "I";
				itGapsText += numGaps;
				itGapsPattern += numGaps;
				//total += numGaps;
				total += 1;
				continue;
			}
			// Count deletions.
			if (isGap(itGapsPattern))
			{
				int numGaps = countGaps(itGapsPattern);
				cigar << numGaps << "D";
				itGapsText += numGaps;
				itGapsPattern += numGaps;
				//total += numGaps;
				total += 1;
				numDels += numGaps;
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

		//result.QEnd = unclippedLength(gapsPattern) - length(gapsPattern) - result.QStart;
		result.QEnd = readLen - (length(gapsPattern) - numDels + result.QStart);
		if(result.QEnd > 0) {
			cigar << result.QEnd << "S";
		}


		//std::cout << "1: " << gapsText << std::endl << "2: " << gapsPattern << std::endl;
		//std::cout << clippedBeginPosition(gapsText) << " " << clippedBeginPosition(gapsPattern) << " " << clippedEndPosition(gapsText) << " " << clippedEndPosition(gapsPattern) << std::endl;
		//std::cout << "Length: " << length(gapsText) << " " << length(gapsPattern) << std::endl;
		//std::cout << "UnClippedLength: " << unclippedLength(gapsText) << " " << unclippedLength(gapsPattern) << std::endl;
		//std::cout << cigar.str() << std::endl;
//
//		if(result.QEnd > 0 || result.QStart > 0) {
//			getchar();
//		}

		//if(cigar.str() == "47S83M1D6M47S") {
		//	getchar();
		//}

		result.Identity = match * 1.0f / total;
		result.NM = mismatch;
		strcpy(result.pBuffer1, cigar.str().c_str());


	}

	virtual int BatchAlign(int const mode, int const batchSize, char const * const * const refSeqList, char const * const * const qrySeqList, Align * const results, void * extData) {

		for(int i = 0; i < batchSize; ++i) {

			TGaps gapsText;
			TGaps gapsPattern;

			//Log.Message("==================================================================================");
			//Log.Message("Ref:  %s", refSeqList[i]);
			//Log.Message("Read: %s", qrySeqList[i]);
			assignSource(gapsText, TSequence(refSeqList[i]));
			assignSource(gapsPattern, TSequence(qrySeqList[i]));

			int len = length(gapsPattern);

			float result = 0.0f;
			switch(m_AlignMode) {
				case 0:
					result = (float)localAlignment(gapsText, gapsPattern, scoringScheme, lDiag, uDiag);
					break;
				case 1:
					result = (float)globalAlignment(gapsText, gapsPattern, scoringScheme, alignConfig, lDiag, uDiag);
					break;
			}
			//float result = (float)globalAlignment(gapsText, gapsPattern, scoringScheme, alignConfig, lDiag, uDiag);
			//float result = (float)localAlignment(gapsText, gapsPattern, scoringScheme, lDiag, uDiag);



			//Log.Message("Score: %f", result);

			//std::cerr << "Score: " << result << std::endl;
			//std::cerr << "The resulting alignment is\n"
			//<< convertToCIGAR(gapsText, gapsPattern) << std::endl;
			//std::cout << "1: " << gapsText << std::endl << "2: " << gapsPattern << std::endl;
			convertToCIGAR(gapsText, gapsPattern, results[i], len);


		}
		return batchSize;
	}

private:
	bool const m_AlignMode;
	int const lDiag;
	int const uDiag;
	Score<float, Simple> scoringScheme;
	seqan::AlignConfig<true, true, false, true> alignConfig;

};

#endif /* ENDTOENDAFFINE_H_ */
