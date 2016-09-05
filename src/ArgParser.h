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

#ifndef __ARGPARSER_H__
#define __ARGPARSER_H__

#include "IConfig.h"

#include <string>

class ArgParser : public IConfig
{
private:

	std::string outDefault;
	std::string noneDefault;

	void ParseArguments(int argc, char const * * argv);

	char * fromString(std::string str);

public:

	ArgParser(int argc, char * argv[]);

	virtual ~ArgParser();

};

#endif //__ARGPARSER_H__
