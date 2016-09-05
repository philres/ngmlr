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

#ifndef __BITPALALIGNER_H__
#define __BITPALALIGNER_H__

#include "IAlignment.h"
#include "IConfig.h"

class BitpalAligner: public IAlignment {

protected:

	int Score(char const * ref, char const * read);

public:
	BitpalAligner();
	~BitpalAligner();

//	void setCorridorSize(int corr);

	int GetScoreBatchSize() const {
		return 8192;
	}
	int GetAlignBatchSize() const {
		return 8192;
	}

	int BatchScore(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList, float * const results,
			void * extData);

	int BatchAlign(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList, Align * const results,
			void * extData);

	int SingleAlign(int const mode, int const corridor,
			char const * const refSeq, char const * const qrySeq,
			Align & result, void * extData);

};

#endif//__BITPALALIGNER_H__
