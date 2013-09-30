#ifndef ___BAMWRITER_H___
#define ___BAMWRITER_H___

#include "GenericReadWriter.h"

#ifdef _BAM

#include "FileWriterBam.h"

class BAMWriter: public GenericReadWriter {
public:
//	BAMWriter(char const * const filename) :
//			GenericReadWriter(filename), file(filename) {
	BAMWriter(FileWriterBam * pWriter, char const * const pFile) :
			GenericReadWriter(), writer(pWriter), file(pFile) {
		NGMInitMutex(&m_OutputMutex);
		bufferIndex = 0;
		//Log.Error("BAM output not supported at the moment!");
		//Fatal();
	}

	~BAMWriter() {
		if (bufferIndex > 0) {
			writer->SaveAlignment(buffer, bufferIndex);
			bufferIndex = 0;
		}
	}

protected:
	virtual void DoWriteProlog();
	virtual void DoWriteRead(MappedRead const * const read, int const scoreId);
	virtual void DoWritePair(MappedRead const * const read1, int const scoreId1, MappedRead const * const read2, int const scoreId2);
	virtual void DoWriteReadGeneric(MappedRead const * const read, int const scoreId, int const pRef, int const pLoc, int const pDist, int const mappingQlty, int flags);
	virtual void DoWriteUnmappedReadGeneric(MappedRead const * const read, int const refId, char const pRefName, int const loc, int const pLoc, int const pDist, int const mappingQlty, int flags);
	virtual void DoWriteUnmappedRead(MappedRead const * const read, int flags = 0x4);
	virtual void DoWriteEpilog();

private:
	char const * const file;
	void translate_flag(BamTools::BamAlignment * al, int flags);

	NGMMutex m_OutputMutex;

	FileWriterBam * writer;

	BamTools::BamAlignment * buffer[10000];

	int bufferIndex;

};
#endif

#endif //___BAMWRITER_H___
