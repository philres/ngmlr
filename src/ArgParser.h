#ifndef __ARGPARSER_H__
#define __ARGPARSER_H__

#include "IConfig.h"


class ArgParser : public IConfig
{
private:

	void ParseArguments(int argc, char const * * argv);

public:

	ArgParser(int argc, char * argv[]);

	virtual ~ArgParser();

};

#endif //__ARGPARSER_H__
