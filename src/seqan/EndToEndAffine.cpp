/*
 * EndToEndAffine.cpp
 *
 *  Created on: Oct 22, 2013
 *      Author: philipp_
 */

#include "EndToEndAffine.h"

int EndToEndAffine::SingleScore(int const mode, int const corridor,
		char const * const refSeq, char const * const qrySeq, float & result,
		void * extData) {

	switch (mode) {
	case 0: {
		TGaps gapsText;
		TGaps gapsPattern;
		assignSource(gapsText, TSequence(refSeq));
		assignSource(gapsPattern, TSequence(qrySeq));
		result = (float) localAlignment(gapsText, gapsPattern, scoringScheme,
				lDiag, uDiag);
	}
		break;
	case 1:
		result = (float) globalAlignmentScore(TSequence(refSeq),
				TSequence(qrySeq), scoringScheme, alignConfig, lDiag, corridor);
		break;
	case 10: {
		TGaps gapsText;
		TGaps gapsPattern;
		assignSource(gapsText, TSequence(refSeq));
		assignSource(gapsPattern, TSequence(qrySeq));
		result = (float) localAlignment(gapsText, gapsPattern, scoringScheme);
	}
		break;
	case 11:
		result = (float) globalAlignmentScore(TSequence(refSeq),
				TSequence(qrySeq), scoringScheme, alignConfig);
		break;
	}
	return 1;
}

int EndToEndAffine::BatchScore(int const mode, int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList, float * const results,
		void * extData) {
	for (int i = 0; i < batchSize; ++i) {
		switch (m_AlignMode) {
		case 0: {
			TGaps gapsText;
			TGaps gapsPattern;
			assignSource(gapsText, TSequence(refSeqList[i]));
			assignSource(gapsPattern, TSequence(qrySeqList[i]));
			results[i] = (float) localAlignment(gapsText, gapsPattern,
					scoringScheme, lDiag, uDiag);
		}
			break;
		case 1:
			results[i] = (float) globalAlignmentScore(TSequence(refSeqList[i]),
					TSequence(qrySeqList[i]), scoringScheme, alignConfig, lDiag,
					uDiag);
			break;
		}
		//results[i] = (float)globalAlignmentScore(TSequence(refSeqList[i]), TSequence(qrySeqList[i]), seqan::MyersBitVector());
	}
	return batchSize;
}

int EndToEndAffine::SingleAlign(int const mode, int const corridor,
		char const * const refSeq, char const * const qrySeq, Align & result,
		void * extData) {

	TGaps gapsText;
	TGaps gapsPattern;
	assignSource(gapsText, TSequence(refSeq));
	assignSource(gapsPattern, TSequence(qrySeq));

	int len = length(gapsPattern);

	float score = 0.0f;
	switch (m_AlignMode) {
	case 0:
		score = (float) localAlignment(gapsText, gapsPattern, scoringScheme,
				lDiag, corridor);
		break;
	case 1:
		score = (float) globalAlignment(gapsText, gapsPattern, scoringScheme,
				alignConfig, lDiag, corridor);
		break;
	case 10:
		score = (float) localAlignment(gapsText, gapsPattern, scoringScheme);
		break;
	case 11:
		score = (float) globalAlignment(gapsText, gapsPattern, scoringScheme,
				alignConfig);
		break;
	}
	convertToCIGAR(gapsText, gapsPattern, result, len);
	result.Score = score;
	return 1;
}

int EndToEndAffine::BatchAlign(int const mode, int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList, Align * const results,
		void * extData) {
	for (int i = 0; i < batchSize; ++i) {

		TGaps gapsText;
		TGaps gapsPattern;
		assignSource(gapsText, TSequence(refSeqList[i]));
		assignSource(gapsPattern, TSequence(qrySeqList[i]));

		int len = length(gapsPattern);

		float result = 0.0f;
		switch (m_AlignMode) {
		case 0:
			result = (float) localAlignment(gapsText, gapsPattern,
					scoringScheme, lDiag, uDiag);
			break;
		case 1:
			result = (float) globalAlignment(gapsText, gapsPattern,
					scoringScheme, alignConfig, lDiag, uDiag);
			break;
		}
		convertToCIGAR(gapsText, gapsPattern, results[i], len);
	}
	return batchSize;
}

void EndToEndAffine::convertToCIGAR(TGaps gapsText, TGaps gapsPattern,
		Align & result, int const readLen) {

	result.QStart = clippedBeginPosition(gapsPattern);
	result.QEnd = 0;
	result.PositionOffset = 0;

	TGapsIterator itGapsPattern = seqan::begin(gapsPattern);
	TGapsIterator itGapsEnd = seqan::end(gapsPattern);

	int match = 0;
	int mismatch = 0;
	int total = 0;

	int count = 0;
	if (m_AlignMode == 1) {
		// Remove trailing gaps in pattern.
		while (isGap(--itGapsEnd))
			++count;
		seqan::setClippedEndPosition(gapsPattern, length(gapsPattern) - count);
	}

	// Remove leading gaps in pattern.
	if (isGap(itGapsPattern)) {
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
	if (result.QStart > 0) {
		cigar << result.QStart << "S";
	}
	int numDels = 0;

	int * errs = new int[10000];
	memset(errs, 0, sizeof(int) * 10000);

	while (itGapsPattern != itGapsEnd) {
		// Count insertions.
		if (isGap(itGapsText)) {
			int numGaps = countGaps(itGapsText);
			cigar << numGaps << "I";
			itGapsText += numGaps;
			itGapsPattern += numGaps;
			//total += numGaps;
			total += 1;
			continue;
		}
		// Count deletions.
		if (isGap(itGapsPattern)) {
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
		while (!isGap(itGapsText) && !isGap(itGapsPattern)
				&& itGapsPattern != itGapsEnd) {
			if (*itGapsText == *itGapsPattern) {
				match += 1;
			} else {
				mismatch += 1;
				int index = total / 100;
				errs[index] += 1;
			}
			total += 1;
			++numChar;
			++itGapsText;
			++itGapsPattern;
		}
		if (numChar != 0) {
			cigar << numChar << "M";
			numChar = 0;
			continue;
		}
	}

	result.QEnd = readLen - (length(gapsPattern) - numDels + result.QStart);
	if (result.QEnd > 0) {
		cigar << result.QEnd << "S";
	}

//		int maxIndex = total / 100;
//		for(int i = 0; i < maxIndex; ++i) {
//			Log.Message("%d: %d", i, errs[i]);
//		}

	delete[] errs;
	result.Identity = match * 1.0f / total;
//		Log.Message("I: %f", result.Identity);
	result.NM = mismatch;
	strcpy(result.pBuffer1, cigar.str().c_str());

}
