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

#include <iostream>
#include <string>

#include <tclap/CmdLine.h>

using std::cerr;
using std::cout;
using std::endl;

class ArgParseOutput : public TCLAP::StdOutput
{
private:

	std::string usageStr;

	std::string versionStr;

public:

	ArgParseOutput(std::string usage, std::string version) {
		usageStr = usage;
		versionStr = version;
	}

	virtual ~ArgParseOutput() {

	}

	virtual void failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e) {
		cerr << "Error:" << endl;
		cerr << "         " << e.error() << endl;
		cerr << endl;
		cerr << "Short usage:" << endl;
		cerr << "        ngmlr [-t <threads>] -r <reference> -q <reads> [-o <output>]" << endl;
		cerr << endl;
		cerr << "For complete USAGE and HELP type:" << endl;
		cerr << "    ngmlr --help" << endl;
		cerr << endl;
		exit(1);
	}

	virtual void usage(TCLAP::CmdLineInterface& c) {
		cerr << usageStr << std::endl;
	}

	virtual void version(TCLAP::CmdLineInterface& c) {

	}
};
