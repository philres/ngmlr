/*
 * OutputReadBuffer.h
 *
 *  Created on: May 3, 2014
 *      Author: philipp_
 */

#ifndef OUTPUTREADBUFFER_H_
#define OUTPUTREADBUFFER_H_

#include <list>

#include "NGMThreads.h"
#include "MappedRead.h"
#include "GenericReadWriter.h"

using std::list;

class OutputReadBuffer {

	static OutputReadBuffer * pInstance;

	NGMMutex m_OutputMutex;

	list<std::pair<MappedRead *, bool> > outputBuffer;

	int currentReadId;

	static int const maxSize = 100000;

	OutputReadBuffer();

	OutputReadBuffer(const OutputReadBuffer& rs) {
		pInstance = rs.pInstance;
		currentReadId = 0;
	}

	OutputReadBuffer& operator =(const OutputReadBuffer& rs) {
		if (this != &rs) {
			pInstance = rs.pInstance;
		}

		return *this;
	}

	~OutputReadBuffer();

public:

	/*Is used by multiple threads but initialized by main thread (before other threads are started), thus no locking is required*/
	static OutputReadBuffer& getInstance() {
		static OutputReadBuffer theInstance;
		pInstance = &theInstance;

		return *pInstance;
	}

	void addRead(MappedRead * read, bool mapped) {
		NGMLock(&m_OutputMutex);

		int index = 0;
		std::list<std::pair<MappedRead *, bool> >::iterator it = outputBuffer.begin();

		while (it != outputBuffer.end() && read->ReadId > (*it).first->ReadId) {
			it++;
		}
		outputBuffer.insert(it, std::pair<MappedRead *, bool>(read, mapped));



		if (outputBuffer.size() > maxSize) {
			throw "Max buffer size reached.";
		}

//		if(outputBuffer.size() >= 2) {
//			for (it=outputBuffer.begin(); it!=outputBuffer.end(); ++it) {
//				std::cout << (*it).first << " " << (*it).first->ReadId << " " << (*it).second << " -- ";
//			}
//			std::cout << std::endl;
//			getchar();
//		}

		NGMUnlock(&m_OutputMutex);
	}

	std::pair<MappedRead *, bool> getNextRead(GenericReadWriter * writer) {
		std::pair<MappedRead *, bool> pair;
		pair.first = 0;
		pair.second = false;
		NGMLock(&m_OutputMutex);

//Log.Message("%d (%s) %d %d", outputBuffer.front().first->ReadId, outputBuffer.front().first->name, currentReadId, outputBuffer.size());

		while (outputBuffer.size() > 0 && outputBuffer.front().first->ReadId == currentReadId) {
//			Log.Message("%s: %d %d", outputBuffer.front().first->name, outputBuffer.front().first->ReadId, currentReadId);
			pair = outputBuffer.front();
			outputBuffer.pop_front();
			currentReadId += 1;
			writer->WriteRead(pair.first, pair.second);
			NGM.GetReadProvider()->DisposeRead(pair.first);
		}

		NGMUnlock(&m_OutputMutex);
		return pair;
	}

};

#endif /* OUTPUTREADBUFFER_H_ */
