/*
 * SWCPU.h
 *
 *  Created on: Jun 15, 2011
 *      Author: fritz
 */

#ifndef CONVEXALIGN_H_
#define CONVEXALIGN_H_

#define pRef pBuffer1
#define pQry pBuffer2

#include "IAlignment.h"
#include "core/Types.h"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

using std::endl;
using std::cout;
using std::max;

#define CIGAR_STOP 10
#define short_min -16000
#define result_number 4
#define line_end '\0'
//#define ref_position 0
//#define qstart 1
//#define qend 2
//#define alignment_offset 3
//#define param_best_read_index 0
//#define param_best_ref_index 1

#define CIGAR_M 0
#define CIGAR_I 1
#define CIGAR_D 2
#define CIGAR_N 3
#define CIGAR_S 4
#define CIGAR_H 5
#define CIGAR_P 6
#define CIGAR_EQ 7
#define CIGAR_X 8

namespace Convex {

class AlignmentMatrix {

public:

	typedef float Score;

	struct MatrixElement {
		Score score;
		short indelRun;
		char direction;

		MatrixElement() {
			score = 0.0f;
			indelRun = 0;
			direction = CIGAR_STOP;
		}
	};

	struct CorridorLine {
		int offset;
		int length;
		ulong offsetInMatrix;
	};

	AlignmentMatrix(int const width, int const height, CorridorLine * corridor,
			int const corridorHeight) {
		privateMatrixHeight = height;
		privateMatrixWidth = width;

		corridorLines = corridor;

		ulong matrixSize = 0;
		for (int i = 0; i < corridorHeight; ++i) {
			corridorLines[i].offsetInMatrix = matrixSize;
			matrixSize += corridorLines[i].length;
		}
		privateMatrix = new MatrixElement[matrixSize];

	}

	virtual ~AlignmentMatrix() {
		if (privateMatrix != 0) {
			delete[] privateMatrix;
			privateMatrix = 0;
		}
	}

	void printMatrix(char const * refSeq, char const * qrySeg) {
		fprintf(stderr, "     - ");
		for (int x = 0; x < privateMatrixWidth; ++x) {
			fprintf(stderr, "  %c ", refSeq[x]);
		}
		fprintf(stderr, "\n");

		for (int y = -1; y < privateMatrixHeight; ++y) {
			if (y == -1) {
				fprintf(stderr, "-: ");
			} else {
				fprintf(stderr, "%c: ", qrySeg[y]);
			}
			for (int x = -1; x < privateMatrixWidth; ++x) {
				fprintf(stderr, "%*d ", 3, (int) getScore(x, y));
			}

			fprintf(stderr, "\n");
		}

	}

	Score getScore(int const x, int const y) {
		MatrixElement * tmp = getElement(x, y);
		return tmp->score;
//		return (tmp == 0) ? 0.0f : tmp->score;
	}

	int getCorridorOffset(int const y) {
		return corridorLines[y].offset;
	}

	int getCorridorLength(int const y) {
		return corridorLines[y].length;
	}

	int getHeight() {
		return privateMatrixHeight;
	}

	int getWidth() {
		return privateMatrixWidth;
	}

	MatrixElement * getElement(int const x, int const y) {
		if (x < 0 || y < 0 || x > (privateMatrixWidth - 1)
				|| y > (privateMatrixHeight - 1)) {
//			empty.score = 0.0f;
			return &empty;
		}
		if (x < corridorLines[y].offset
				|| x >= (corridorLines[y].offset + corridorLines[y].length)) {
//			empty.score = 0.0f;
			return &empty;
		}

		return privateMatrix + corridorLines[y].offsetInMatrix
				+ (x - corridorLines[y].offset);
	}

private:

	int privateMatrixHeight;

	int privateMatrixWidth;

	MatrixElement * privateMatrix;

	MatrixElement empty;

	CorridorLine * corridorLines;

};

//const char trans[256] = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0,
//		4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4,
//		4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

class ConvexAlign: public IAlignment {
public:
	ConvexAlign(int gpu_id);
	virtual ~ConvexAlign();
	virtual int GetScoreBatchSize() const;
	virtual int GetAlignBatchSize() const;
	virtual int BatchScore(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList, float * const results,
			void * extData);
	virtual int BatchAlign(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList, Align * const results,
			void * extData);

	virtual int SingleAlign(int const mode, int const corridor,
			char const * const refSeq, char const * const qrySeq,
			Align & result, void * extData);
private:

	struct FwdResults {
		int best_ref_index;
		int best_read_index;
		int max_score;
		int qend;
		int qstart;
		int ref_position;
		int alignment_offset;
	};

	//bool cigar;
	//short scores[6][6];
	AlignmentMatrix::Score mat;
	AlignmentMatrix::Score mis;
	AlignmentMatrix::Score gap_open_read;
	AlignmentMatrix::Score gap_open_ref;
	AlignmentMatrix::Score gap_ext;
	AlignmentMatrix::Score gap_ext_min;
	AlignmentMatrix::Score gap_decay;

	AlignmentMatrix * matrix;

//	long long maxAlignMatrixLen;
//	MatrixElement * alignMatrix;
	int * binaryCigar;

	//meta info
	unsigned int batch_size; //effictive thread number that is started per call

	int const maxBinaryCigarLength;

	bool const pacbioDebug;

	int printCigarElement(char const op, int const length, char * cigar);

	int convertCigar(char const * const refSeq, Align & result,
			FwdResults & fwdResults, int const externalQStart,
			int const externalQEnd);

	AlignmentMatrix::Score fwdFillMatrix(char const * const refSeq,
			char const * const qrySeq, FwdResults & fwdResults,
			int readId);

	bool revBacktrack(char const * const refSeq, char const * const qrySeq,
			FwdResults & fwdResults, int readId);

};

#endif /* CONVEXALIGN_H_ */

}  // namespace Convex
