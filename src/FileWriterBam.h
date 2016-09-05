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

#ifndef FILEWRITERBAM_H_
#define FILEWRITERBAM_H_

#include <string.h>

#include "api/BamAlignment.h"
#include "api/BamWriter.h"
#include "api/SamHeader.h"

#include "NGMThreads.h"
#include "ILog.h"
#include "IConfig.h"

static int countBam = 0;

class FileWriterBam {

public:

	NGMMutex m_OutputMutex;

	FileWriterBam(char const * const filename) {
		countBam += 1;
		if (countBam > 1) {
			Log.Error("To many FileWriter instances!");
		}
		NGMInitMutex(&m_OutputMutex);
		writer = new BamTools::BamWriter();
	}

	void SaveAlignment(BamTools::BamAlignment * buffer[], int const n) {
		NGMLock(&m_OutputMutex);

		for (int j = 0; j < n; ++j) {
			if (!writer->SaveAlignment(*buffer[j])) {
				Log.Error("Couldn't write BAM record!");
			}
			delete buffer[j];
			buffer[j] = 0;
		}
		NGMUnlock(&m_OutputMutex);
	}

	bool Open(char const * const filename, BamTools::SamHeader header, BamTools::RefVector refs) {
		if(filename == 0) {
			return writer->Open("stdout", header, refs);
		} else {
			return writer->Open(filename, header, refs);
		}
	}

	~FileWriterBam() {
		NGMLock(&m_OutputMutex);
		if (writer != 0) {
			writer->Close();
			delete writer;
			writer = 0;
		}
		NGMUnlock(&m_OutputMutex);

	}

private:
	BamTools::BamWriter * writer;

};

#endif /* FILEWRITERBAM_H_ */
