/*
 * oclTool.cpp
 *
 *  Created on: Sep 22, 2011
 *      Author: philipp_
 */

#include <iostream>
#include <fstream>
#include <string>

using std::ifstream;
using std::string;
using std::cout;
using std::ofstream;

int main(int argc, char *argv[]) {
	string line;
	ifstream myfile;

	string outputName(argv[argc - 1]);
	if (outputName.find(".cl") == string::npos) {

		FILE * pFile;
		pFile = fopen(outputName.c_str(), "w");

		if (pFile != 0) {
			fprintf(pFile, "char const %s[] = {\n", argv[1]);
			for (int i = 2; i < argc - 1; ++i) {
				myfile.open(argv[i]);
				if (myfile.is_open()) {
					while (myfile.good()) {
						getline(myfile, line);
						for (int j = 0; j < line.length(); ++j) {
							fprintf(pFile, "0x%02x, ", (unsigned char) line[j]);
						}
						fprintf(pFile, "0x%02x, ", '\n');
					}
					myfile.close();
				}
				myfile.close();
			}
			fprintf(pFile, "0x%02x };\n", '\0');
			fclose (pFile);
		} else
			cout << "Unable to open file";
	} else {
		cout << "Stop!" << std::endl;
	}

	return 0;
}
