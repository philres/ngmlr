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

#include "AlignmentMatrixFast.h"

#include <iostream>
#include <stdio.h>

using std::cerr;

namespace Convex {

AlignmentMatrixFast::AlignmentMatrixFast() :
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

void AlignmentMatrixFast::prepare(int const width, int const height,
		CorridorLine * corridor, int const corridorHeight) {
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

	return x > minCorridor && x < maxCorridor;
}

}  // namespace Convex
