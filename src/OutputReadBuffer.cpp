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
