/**
 * Contact: philipp.rescheneder@gmail.com
 */

#include "AlignmentMatrixFast.h"

#include <iostream>
#include <stdio.h>

using std::cerr;

namespace Convex {

AlignmentMatrixFast::AlignmentMatrixFast(ulong const pMaxMatrixSizeMB) :
		STOP(CIGAR_STOP), defaultMatrixSize((ulong) 30000 * (ulong) 9000), maxMatrixSizeMB(pMaxMatrixSizeMB) {
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

bool AlignmentMatrixFast::prepare(int const width, int const height,
		CorridorLine * corridor, int const corridorHeight) {
	privateMatrixHeight = height;
	privateMatrixWidth = width;

	corridorLines = corridor;

	uloc matrixSize = 0;
	int maxCorridorLength = 0;
	for (int i = 0; i < corridorHeight; ++i) {
		corridorLines[i].offsetInMatrix = matrixSize;
		matrixSize += corridorLines[i].length;
		maxCorridorLength = std::max(maxCorridorLength,
				corridorLines[i].length);
	}
	if ((ulong) (matrixSize * sizeof(char) / 1000.0f / 1000.0f) < maxMatrixSizeMB) {
		if (matrixSize > privateMatrixSize) {
//			fprintf(stderr, "Reallocating matrix for alignment (size %llu)\n\n", (long long) (matrixSize * sizeof(char) / 1000.0f / 1000.0f));
			delete[] directionMatrix;
			directionMatrix = 0;
			directionMatrix = new char[matrixSize];
			privateMatrixSize = matrixSize;
		}
		currentLine = new MatrixElement[maxCorridorLength];
		lastLine = new MatrixElement[maxCorridorLength];
	} else {
		fprintf(stderr, "Warning: Couldn't allocate alignment matrix. Required memory (%llu) > max matrix size (%lu)\n\n", (long long) (matrixSize * sizeof(char) / 1000.0f / 1000.0f), maxMatrixSizeMB);
		return false;
	}
	return true;
}

void AlignmentMatrixFast::clean() {
	if (currentLine != 0) {
		delete[] currentLine;
		currentLine = 0;
	}
	if (lastLine != 0) {
		delete[] lastLine;
		lastLine = 0;
	}
	if (privateMatrixSize > defaultMatrixSize) {
		delete[] directionMatrix;
		privateMatrixSize = defaultMatrixSize;
		directionMatrix = new char[privateMatrixSize];
	}
}

AlignmentMatrixFast::~AlignmentMatrixFast() {
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

void AlignmentMatrixFast::printMatrix(char const * refSeq, char const * qrySeg) {
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

AlignmentMatrixFast::Score AlignmentMatrixFast::getScore(int const x, int const y) {
	return getElement(x, y)->score;
}

int AlignmentMatrixFast::getCorridorOffset(int const y) const {
	return corridorLines[y].offset;
}

int AlignmentMatrixFast::getCorridorLength(int const y) const {
	return corridorLines[y].length;
}

int AlignmentMatrixFast::getHeight() const {
	return privateMatrixHeight;
}

int AlignmentMatrixFast::getWidth() const {
	return privateMatrixWidth;
}

AlignmentMatrixFast::MatrixElement * AlignmentMatrixFast::getElement(int const x, int const y) {
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

AlignmentMatrixFast::MatrixElement * AlignmentMatrixFast::getElementEdit(int const x, int const y) {
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

char * AlignmentMatrixFast::getDirection(int const x, int const y) {
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

bool AlignmentMatrixFast::prepareLine(int const y) {
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

bool AlignmentMatrixFast::validPath(int const x, int const y) {
	int width = corridorLines[y].length;
	int minCorridor = corridorLines[y].offset + 0.1f * width;
	int maxCorridor = minCorridor + width - 0.1f * width;
//	fprintf(stderr, "offest: %d, length: %d\n", corridorLines[y].offset, width);
//	fprintf(stderr, "%d < %d < %d\n", minCorridor, x, maxCorridor);
	return x > minCorridor && x < maxCorridor;
}

}  // namespace Convex
