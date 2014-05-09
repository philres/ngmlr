/*
 * OutputReadBuffer.cpp
 *
 *  Created on: May 3, 2014
 *      Author: philipp_
 */

#include "OutputReadBuffer.h"

OutputReadBuffer * OutputReadBuffer::pInstance = 0;

OutputReadBuffer::OutputReadBuffer() {
	NGMInitMutex(&m_OutputMutex);
	currentReadId = 0;

}

OutputReadBuffer::~OutputReadBuffer() {
	if(outputBuffer.size() > 0) {
		throw "Elements left in buffer!";
	}
}

//#include <iostream>
//
////g++ ../../src/OutputReadBuffer.cpp ../../src/core/unix_threads.cpp ../../src/MappedRead.cpp -I ../../src/core/ -I ../../include/ -o buffer -pthread
//
//void write(OutputReadBuffer & buffer) {
//	std::pair<MappedRead *, bool> pair;
//	while ((pair = buffer.getNextRead()).first != 0) {
//		std::cout << "Printing read: " << pair.first->ReadId << std::endl;
//	}
//	std::cout << "---------------" << std::endl;
//}
//
//int main() {
//	OutputReadBuffer & buffer = OutputReadBuffer::getInstance();
//
//	buffer.addRead(new MappedRead(4, 100), false);
//	write(buffer);
//	buffer.addRead(new MappedRead(1, 100), false);
//	write(buffer);
//	buffer.addRead(new MappedRead(0, 100), false);
//	write(buffer);
//	buffer.addRead(new MappedRead(2, 100), false);
//	write(buffer);
//
//
//	buffer.addRead(new MappedRead(3, 100), false);
//	write(buffer);
//	buffer.addRead(new MappedRead(5, 100), false);
//	write(buffer);
//	buffer.addRead(new MappedRead(6, 100), false);
//	write(buffer);
//
//
//
//	return 0;
//}
