/*
 * FileWriterBam.h
 *
 *  Created on: Mar 2, 2013
 *      Author: philipp
 */

#ifndef FILEWRITERBAM_H_
#define FILEWRITERBAM_H_

#include "api/BamAlignment.h"
#include "api/BamWriter.h"
#include "api/SamHeader.h"

#include "NGMThreads.h"
#include "ILog.h"
#include "Config.h"

static int countBam = 0;

class FileWriterBam {

public:

	NGMMutex m_OutputMutex;

	FileWriterBam(char const * const filename) {
		countBam += 1;
		if (countBam > 1) {
			Log.Error("To many FileWriter instances!");
			Fatal();
		}
		NGMInitMutex(&m_OutputMutex);
		writer = new BamTools::BamWriter();
		//i = 0;
	}

	void SaveAlignment(BamTools::BamAlignment * buffer[], int const n) {
		NGMLock(&m_OutputMutex);
//		if (i < 10000) {
//			buffer[i++] = al;
//		} else {
			for (int j = 0; j < n; ++j) {
				if (!writer->SaveAlignment(*buffer[j])) {
					Log.Error("Couldn't write BAM record!");
					Fatal();
				}
				delete buffer[j]; buffer[j] = 0;
			}
//			i = 0;
//		}
		NGMUnlock(&m_OutputMutex);
	}

	bool Open(char const * const filename, BamTools::SamHeader header, BamTools::RefVector refs) {
		return writer->Open(filename, header, refs);
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

//	BamTools::BamAlignment buffer[10000];
//
//	int i;

};

#endif /* FILEWRITERBAM_H_ */
