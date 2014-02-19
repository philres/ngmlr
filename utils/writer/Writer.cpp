/*
 * NgmWriter.cpp
 *
 *  Created on: May 4, 2011
 *      Author: philipp_
 */

#include "Writer.h"

#include <string.h>
#include <stdarg.h>
#include <ctype.h>

#include "Log.h"

#undef module_name
#define module_name "WRITER"

int Writer::Print(const char *format, ...) {
	va_list arg;
	int done;

	va_start(arg, format);
	done = vsprintf(writeBuffer + bufferPosition, format, arg);
	bufferPosition += done;
	va_end(arg);
	return done;
}

void Writer::Flush(bool last) {
	if (bufferPosition > BUFFER_LIMIT || last) {
		if(fwrite(writeBuffer, sizeof(char), bufferPosition, m_Output) != bufferPosition) {
			throw "Couldn't write to file";
		}
		bufferPosition = 0;
		fflush(m_Output);
	}
}

Writer::Writer(char const * const fileName) {
	if (!(m_Output = fopen(fileName, "wb"))) {
		Log.Error("Unable to open output file %s", fileName);
		Fatal();
	}
	writeBuffer = new char[BUFFER_SIZE];
	bufferPosition = 0;
	Log.Message("Opening %s for writing.", fileName);
}

void Writer::toUpperCase(char * sequence, int sequenceLength) {
	int i = 0;
	while (sequence[i] != '\0') {
		sequence[i] = ::toupper(sequence[i]);
		i += 1;
	}
}

Writer::~Writer() {
	Flush(true);
	fclose(m_Output);
}
