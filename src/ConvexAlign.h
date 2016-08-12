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

	AlignmentMatrix() :
			STOP(CIGAR_STOP), defaultMatrixSize((ulong) 30000 * (ulong) 9000) {
		privateMatrixHeight = 0;
		privateMatrixWidth = 0;
		corridorLines = 0;

		privateMatrixSize = defaultMatrixSize;
		directionMatrix = new char[privateMatrixSize];
		currentLine = 0;
		lastLine = 0;

		lastY = 0;
		currentY = 0;

	}

	void prepare(int const width, int const height, CorridorLine * corridor,
			int const corridorHeight) {
		privateMatrixHeight = height;
		privateMatrixWidth = width;

		corridorLines = corridor;

		ulong matrixSize = 0;
		int maxCorridorLength = 0;
		for (int i = 0; i < corridorHeight; ++i) {
			corridorLines[i].offsetInMatrix = matrixSize;
			matrixSize += corridorLines[i].length;
			maxCorridorLength = std::max(maxCorridorLength,
					corridorLines[i].length);
		}
		if (matrixSize > privateMatrixSize) {
			//fprintf(stderr, "Reallocating matrix for alignment\n");
			delete[] directionMatrix;
			directionMatrix = 0;
			directionMatrix = new char[matrixSize];
			privateMatrixSize = matrixSize;
		}

		currentLine = new MatrixElement[maxCorridorLength];
		lastLine = new MatrixElement[maxCorridorLength];

	}

	void clean() {
		if (currentLine != 0) {
			delete[] currentLine;
			currentLine = 0;
		}
		if (lastLine != 0) {
			delete[] lastLine;
			lastLine = 0;
		}
		if(privateMatrixSize > defaultMatrixSize) {
			delete[] directionMatrix;
			privateMatrixSize = defaultMatrixSize;
			directionMatrix = new char[privateMatrixSize];
		}
	}

	virtual ~AlignmentMatrix() {
		if (directionMatrix != 0) {
			delete[] directionMatrix;
			directionMatrix = 0;
		}
		if (currentLine != 0) {
			delete[] currentLine;
			currentLine = 0;
		}
		if (lastLine != 0) {
			delete[] lastLine;
			lastLine = 0;
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
				fprintf(stderr, "%*d ", 3, (int) *getDirection(x, y));
			}

			fprintf(stderr, "\n");
		}

	}

	Score getScore(int const x, int const y) {
		return getElement(x, y)->score;
	}

	int getCorridorOffset(int const y) const {
		return corridorLines[y].offset;
	}

	int getCorridorLength(int const y) const {
		return corridorLines[y].length;
	}

	int getHeight() const {
		return privateMatrixHeight;
	}

	int getWidth() const {
		return privateMatrixWidth;
	}

	MatrixElement * getElement(int const x, int const y) {
		if (y < 0 || x < 0) {
			return &empty;
		}
		if (y != currentY && y != lastY) {
			fprintf(stderr, "Index not found\n");
			throw "Index not found";
		} else if (y == currentY) {
			if (x < currentCorridor.offset
					|| x >= (currentCorridor.offset + currentCorridor.length)) {
				return &empty;
			}
			return currentLine + (x - currentCorridor.offset);
		} else { // y == lastY
			if (x < lastCorridor.offset
					|| x >= (lastCorridor.offset + lastCorridor.length)) {
				return &empty;
			}
			return lastLine + (x - lastCorridor.offset);
		}
		return 0;
	}

	MatrixElement * getElementEdit(int const x, int const y) {
		if (y < 0 || x < 0) {
			fprintf(stderr, "Element not found in alignment matrix.\n");
			throw "";
		}
		if (y != currentY && y != lastY) {
			fprintf(stderr, "Index not found\n");
			throw "Index not found";
		} else if (y == currentY) {
			if (x < currentCorridor.offset
					|| x >= (currentCorridor.offset + currentCorridor.length)) {
				fprintf(stderr, "Element not found in alignment matrix.\n");
				throw "";
			}
			return currentLine + (x - currentCorridor.offset);
		} else { // y == lastY
			if (x < lastCorridor.offset
					|| x >= (lastCorridor.offset + lastCorridor.length)) {
				fprintf(stderr, "Element not found in alignment matrix.\n");
				throw "";
			}
			return lastLine + (x - lastCorridor.offset);
		}

		return 0;
	}

	char * getDirection(int const x, int const y) {
		if (y < 0 || y > (privateMatrixHeight - 1) || x < 0) {
			return &STOP;
		} else {
			CorridorLine line = corridorLines[y];
			if (x < line.offset || x >= (line.offset + line.length)) {
				return &STOP;
			}
			return directionMatrix + line.offsetInMatrix + (x - line.offset);
		}
	}

	bool prepareLine(int const y) {
		if (y < 0 || y > getHeight()) {
			return false;
		}
		MatrixElement * tmp = lastLine;

		lastLine = currentLine;
		lastY = currentY;
		lastCorridor = currentCorridor;

		currentLine = tmp;
		currentY = y;
		currentCorridor = corridorLines[currentY];
		return true;
	}

	bool validPath(int const x, int const y) {
		int width = corridorLines[y].length;
		int minCorridor = corridorLines[y].offset + 0.1f * width;
		int maxCorridor = minCorridor + width - 0.1f * width;

		return x > minCorridor && x < maxCorridor;
	}

private:

	ulong const defaultMatrixSize;

	int privateMatrixHeight;

	int privateMatrixWidth;

	ulong privateMatrixSize;

	MatrixElement * currentLine;

	CorridorLine currentCorridor;

	int currentY;

	MatrixElement * lastLine;

	CorridorLine lastCorridor;

	int lastY;

	char * directionMatrix;

	MatrixElement empty;

	char STOP;

	CorridorLine * corridorLines;

};

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

	virtual int SingleAlign(int const mode, CorridorLine * corridor,
			int const corridorHeight, char const * const refSeq,
			char const * const qrySeq, Align & result, int const externalQStart,
			int const externalQEnd, void * extData);

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

	AlignmentMatrix::Score mat;
	AlignmentMatrix::Score mis;
	AlignmentMatrix::Score gap_open_read;
	AlignmentMatrix::Score gap_open_ref;
	AlignmentMatrix::Score gap_ext;
	AlignmentMatrix::Score gap_ext_min;
	AlignmentMatrix::Score gap_decay;

	AlignmentMatrix * matrix;

	int * binaryCigar;

	//meta info
	unsigned int batch_size; //effictive thread number that is started per call

	int const maxBinaryCigarLength;

	bool const pacbioDebug;

	bool const stdoutPrintAlignCorridor;

	int printCigarElement(char const op, int const length, char * cigar);

	int convertCigar(char const * const refSeq, Align & result,
			FwdResults & fwdResults, int const externalQStart,
			int const externalQEnd);

	AlignmentMatrix::Score fwdFillMatrix(char const * const refSeq,
			char const * const qrySeq, FwdResults & fwdResults, int readId);

	bool revBacktrack(char const * const refSeq, char const * const qrySeq,
			FwdResults & fwdResults, int readId);

};

#endif /* CONVEXALIGN_H_ */

}  // namespace Convex
