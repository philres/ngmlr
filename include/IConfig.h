#ifndef __ICONFIG_H__
#define __ICONFIG_H__

class IConfig
{
public:
	virtual char const * GetString(char const * const name) const = 0;
	virtual int GetInt(char const * const name) const = 0;
	virtual int GetInt(char const * const name, int min, int max) const = 0;
	virtual int GetParameter(char const * const name) const = 0;
	virtual float GetFloat(char const * const name) const = 0;
	virtual float GetFloat(char const * const name, float min, float max) const = 0;

	virtual int GetIntArray(char const * const name, int * pData, int len) const = 0;
	virtual int GetFloatArray(char const * const name, float * pData, int len) const = 0;
	virtual int GetDoubleArray(char const * const name, double * pData, int len) const = 0;

	virtual bool Exists(char const * const name) const = 0;
	virtual bool HasArray(char const * const name) const = 0;

	virtual ~IConfig() {};
};

typedef void (*pfSetConfig)(IConfig const *);

extern IConfig* _config;
#define Config (*_config)

#endif
