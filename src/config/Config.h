#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "IConfig.h"

#include <map>
#include <string>

class _Config : public IConfig
{
public:
	char const * GetString(char const * const name) const;
	int GetInt(char const * const name) const;
	int GetInt(char const * const name, int min, int max) const;
	int GetParameter(char const * const name) const;
	float GetFloat(char const * const name) const;
	float GetFloat(char const * const name, float min, float max) const;

	int GetIntArray(char const * const name, int * pData, int len) const;
	int GetFloatArray(char const * const name, float * pData, int len) const;
	int GetDoubleArray(char const * const name, double * pData, int len) const;

	bool Exists(char const * const name) const;
	bool HasArray(char const * const name) const;

	void Override(char const * const name, char const * const value);
	void Override(char const * const name, int const value);
	void Override(char const * const name, float const value);

	void Default(char const * const name, char const * const value);
	void Default(char const * const name, int const value);
	void Default(char const * const name, float const value);

	_Config(int argc, char * argv[], bool praseArgs = true);
	virtual ~_Config();
private:
	bool initialized;
	bool sealed;

	typedef std::map<std::string, std::string> TConfigMap;
	TConfigMap * config_map;
	TConfigMap * config_arrays;

	bool InternalExists(std::string name) const;
	std::string InternalGet(std::string name, char const * * arr_data) const;
	std::string InternalGet(std::string name) const;
	void InternalAdd(std::string name, std::string value, std::string arr_data, bool override = false);
	void ParseFile(char const * const filename);
	void ParseArguments(int argc, char * argv[]);

	template <typename T> friend int ParseMatrix(_Config const * config, char const * const name, T * data, int len, T f(char const * const));
};

#endif
