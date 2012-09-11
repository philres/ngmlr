#ifndef __GENERICREADWRITER_H__
#define __GENERICREADWRITER_H__

#include <stdio.h>
#include <stdarg.h>

#include "ILog.h"
#include "Config.h"

#include "MappedRead.h"

class GenericReadWriter {
public:
	GenericReadWriter(char const * const filename) {
		if (!(m_Output = fopen(filename, "wb"))) {
			Log.Error("Unable to open output file %s", filename);
		}
		writeBuffer = new char[BUFFER_SIZE];
		bufferPosition = 0;

		identity = 0;
		if (Config.Exists("identity_tresh")) {
			identity = Config.GetFloat("identity_tresh", 0, 100);
			identity = identity / 100.0f;
		}
	}
	virtual ~GenericReadWriter() {
		if (m_Output) {
			Flush(true);
			fclose(m_Output);
		}
	}
protected:

	virtual void DoWriteProlog() = 0;
	virtual void DoWriteRead(MappedRead const * const read) = 0;
	virtual void DoWritePair(MappedRead const * const read1, MappedRead const * const read2) = 0;
	virtual void DoWriteUnmappedRead(MappedRead const * const read) = 0;
	virtual void DoWriteEpilog() = 0;

	int Print(const char *format, ...) {
		va_list arg;
		int done;

		va_start(arg, format);
		done = vsprintf(writeBuffer + bufferPosition, format, arg);
		bufferPosition += done;
		va_end(arg);
		return done;
	}

	void Flush(bool last = false) {
		if (bufferPosition > BUFFER_LIMIT || last) {
			fwrite(writeBuffer, sizeof(char), bufferPosition, m_Output);
			bufferPosition = 0;
			fflush(m_Output);
		}
	}

	float identity;

private:
	FILE * m_Output;

	static int const BUFFER_SIZE = 17000000;
	static int const BUFFER_LIMIT = 16000000;

	char * writeBuffer;
	int bufferPosition;

public:
	void WriteProlog() {
		if (m_Output)
			DoWriteProlog();
	}

	void WriteRead(MappedRead const * const read, bool mapped = true) {
		if (m_Output) {

			static float const minIdentity = Config.GetFloat("min_identity", 0.0f, 1.0f);
			static int const minResidues = Config.GetInt("min_residues", 0, 1000);

			mapped = mapped && (read->Identity >= minIdentity);
			mapped = mapped && ((read->length - read->QStart - read->QEnd) >= minResidues);

			if (mapped) {
				DoWriteRead(read);
			} else {
				DoWriteUnmappedRead(read);
			}
			Flush();
		}
	}

	void WritePair(MappedRead * const read1, MappedRead * const read2) {
		if (m_Output) {
			static float const minIdentity = Config.GetFloat("min_identity", 0.0f, 1.0f);
			static int const minResidues = Config.GetInt("min_residues", 0, 1000);

			bool mapped1 = read1->hasCandidates() && (read1->Identity >= minIdentity)
					&& ((read1->length - read1->QStart - read1->QEnd) >= minResidues);
			bool mapped2 = read2->hasCandidates() && (read2->Identity >= minIdentity)
					&& ((read2->length - read2->QStart - read2->QEnd) >= minResidues);



			if (!mapped1) {
				read1->clearScores(false);
			}
			if (!mapped2) {
				read2->clearScores(false);
			}
			//Log.Verbose("Output paired 1: hC %d, R: %d %d, I: %f %f", read1->hasCandidates(), ((read1->length - read1->QStart - read1->QEnd)), minResidues, read1->Identity >= minIdentity, minIdentity);
			//Log.Verbose("%d %d %d", read1->length, read1->QStart, read1->QEnd);
			//Log.Verbose("Output paired 2: hC %d, R: %d %d, I: %f %f", read2->hasCandidates(), ((read2->length - read2->QStart - read2->QEnd)), minResidues, read2->Identity >= minIdentity, minIdentity);
			//Log.Verbose("%d %d %d", read2->length, read2->QStart, read2->QEnd);

			DoWritePair(read1, read2);

			Flush();
		}
	}

//	void WriteRead(std::list<MappedRead const * const> reads) {
//		if (m_Output) {
//			for(int i = 0; i < reads.size(); ++i) {
//				WriteRead(reads[i]);
//			}
//		}
//	}
	void WriteEpilog() {
		if (m_Output)
			DoWriteEpilog();
	}
};

#endif
