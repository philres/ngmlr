/**
 * Contact: philipp.rescheneder@gmail.com
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
