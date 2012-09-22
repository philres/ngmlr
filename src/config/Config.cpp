#include "Config.h"

#include "PlatformSpecifics.h"
#include "Log.h"
#include "NGM.h"
#include "Debug.h"
#include "Options.h"

#include <map>
#include <string>
#include <algorithm>
#include <sstream>

#include <memory.h>
#include <ctype.h>
#include <unistd.h>
#include <getopt.h>
#include <stdio.h>     /* for printf */
#include <stdlib.h>    /* for exit */

//static const struct option long_options[] = {
//		{ "config", required_argument, 0, 'c' },
//		{ "output", no_argument, 0, '0' },
//		{ "CCC", required_argument, 0, 'c' },
//		{ "DDD", no_argument, 0, 'd' },
//		{ "EEE", required_argument, 0, 0 },
//		{ "FFF", required_argument, 0, 0 }, 0 };

extern char *optarg;
extern int optind, opterr, optopt;

#undef module_name
#define module_name "CONFIG"

typedef std::map<std::string, std::string> SwitchMap;

const SwitchMap::value_type rawData[] = { SwitchMap::value_type("c", "config"), SwitchMap::value_type("o", "output"), };
const int numElems = sizeof rawData / sizeof rawData[0];
SwitchMap switchMap(rawData, rawData + numElems);

std::string add_timestamp(std::string); // from logging.cpp

// Config is stored as a string to string map, with
// values converted just in time for the user to
// provide type information rather then trying to find
// out ourselves

bool _Config::InternalExists(std::string name) const {
	std::transform(name.begin(), name.end(), name.begin(), ::tolower);
	return config_map->count(name) == 1;
}

std::string _Config::InternalGet(std::string name, char const * * arr_data) const {
	std::transform(name.begin(), name.end(), name.begin(), ::tolower);
	if (!InternalExists(name)) {
		Log.Warning("Tried to access unknown config value \"%s\"", name.c_str());
		return "";
	} else {
		if (arr_data != 0) {
			if (config_arrays->count(name) == 0) {
				Log.Error("No array declared for %s", name.c_str());
			} else {
				*arr_data = (*config_arrays)[name].c_str();
			}
		}
		return (*config_map)[name];
	}
}
std::string _Config::InternalGet(std::string name) const {
	return InternalGet(name, 0);
}

int replace(int i) {
	return (i == '-') ? '_' : i;
}

void _Config::InternalAdd(std::string name, std::string value, std::string arr_data = std::string(), bool override) {
	std::transform(name.begin(), name.end(), name.begin(), ::tolower);
	std::transform(name.begin(), name.end(), name.begin(), replace);

	bool exists = InternalExists(name);
	if (exists && !override) {
		Log.Warning("Ignoring redefinition of %s (current value = %s, tried to override with %s)", name.c_str(), InternalGet(name).c_str(), value.c_str());
	} else {
		if (name == "logfile" || name == "output")
			value = add_timestamp(value);

		if (name == "qry_end")
			Log.Warning("Parameter qry_end has been replaced by qry_count");

		if (name != "cmdline") {
			Log.Verbose("%s <%s : %s>%s", (exists) ? "Overriding" : "Adding", name.c_str(), value.c_str(), (arr_data.length() != 0) ? "[Array]" : "");
		}
		config_map->insert(std::pair<std::string, std::string>(name, value));
		if (arr_data.length() != 0)
			config_arrays->insert(std::pair<std::string, std::string>(name, arr_data));
	}
}

extern char const * const defaultConfigText;

void CreateDefaultConfig() {
	Log.Message("Writing default config");

	FILE * f = fopen("config.default", "w");
	if (f == 0) {
		Log.Error("Unable to create file");
		Fatal();
	}
	fprintf(f, "%s", defaultConfigText);
	fclose(f);
	exit(0);
}

void ParseArguments(int argc, char * argv[]);

/*_Config const & _Config::Instance()
 {
 if (pInstance == 0)
 {
 pInstance = new _Config();
 }

 return *pInstance;
 }*/

char const * _Config::GetString(char const * const name) const {
	return InternalGet(name).c_str();
}

int GetConversionFactor(std::string s) {
	char c = s[s.length() - 1];
	if (c == 'k')
		return 1000;
	else if (c == 'm')
		return 1000 * 1000;
	else if (c == 'g')
		return 1000 * 1000 * 1000;
	else
		return 1;
}

int _Config::GetInt(char const * const name) const {
	std::string const str = InternalGet(name);
	if (!str.empty())
		return atoi(str.c_str()) * GetConversionFactor(str);
	else
		return 0;
}

int _Config::GetParameter(char const * const name) const {
	return GetInt(name);
}

float _Config::GetFloat(char const * const name) const {
	std::string const str = InternalGet(name);
	if (!str.empty())
		return (float) atof(str.c_str()) * (float) GetConversionFactor(str);
	else
		return 0.0f;
}

int _Config::GetInt(char const * const name, int min, int max) const {
	int i = GetInt(name);
	if ((i < min) || ((max >= min) && (i > max))) {
		if (max >= min)
			Log.Error("Value %s : %i out of range [%i, %i] - defaulting to %i", name, i, min, max, min);
		else
			Log.Error("Value %s : %i below minimum of %i - defaulting to %i", name, i, min, min);

		i = min;
	}
	return i;
}
float _Config::GetFloat(char const * const name, float min, float max) const {
	float f = GetFloat(name);
	if ((f < min) || ((max >= min) && (f > max))) {
		Log.Error("Value %s : %f out of range [%f, %f]", name, f, min, max);
		Fatal();
	}
	return f;
}

void _Config::Override(char const * const name, char const * const value) {
	InternalAdd(name, value, "", true);
}
void _Config::Override(char const * const name, int const value) {
	char buffer[32];
	sprintf(buffer, "%i", value);
	Override(name, buffer);
}

void _Config::Override(char const * const name, float const value) {
	char buffer[32];
	sprintf(buffer, "%f", value);
	Override(name, buffer);
}

void _Config::Default(char const * const name, char const * const value) {
	if (!InternalExists(name)) {
		Log.Verbose("Parameter %s not set. Using default value (%s).", name, value);
		InternalAdd(name, value, "", false);
	}
}
void _Config::Default(char const * const name, int const value) {
	if (!InternalExists(name)) {
		Log.Verbose("Parameter %s not set. Using default value (%d).", name, value);
		char buffer[32];
		sprintf(buffer, "%i", value);
		InternalAdd(name, buffer, "", false);
	}
}

void _Config::Default(char const * const name, float const value) {
	if (!InternalExists(name)) {
		Log.Verbose("Parameter %s not set. Using default value (%f).", name, value);
		char buffer[32];
		sprintf(buffer, "%f", value);
		InternalAdd(name, buffer, "", false);
	}
}

// Helper functions for config file parsing
inline bool IsWhitespace(char const * const str) {
	return (*str == ' ' || *str == '\n' || *str == '\t' || *str == '\r');
}
inline void SkipNonWhitespace(char const * & str) {
	while (*str != 0 && !IsWhitespace(str))
		++str;
}
inline void SkipWhitespace(char const * & str) {
	while (*str != 0 && IsWhitespace(str))
		++str;
}
inline void SkipLine(char const * & str, bool stopForArray) {
	while (*str != 0 && *str != '\n' && *str != '\r' && !(stopForArray && *str == '{')) {
		++str;
	}
	if (*str != '{')
		++str;
}

// Converts the stored matrix into actual numbers by the given conversion function
// and writes them into the data array
template<typename T> int ParseMatrix(_Config const * config, char const * const name, T * data, int len, T f(char const * const)) {
	char const * matrix_str = 0;

	char const * arr_size_str = config->InternalGet(name, &matrix_str).c_str();
	if (arr_size_str == 0) {
		Log.Error("Invalid matrix \"%s\"", name);
		Fatal();
	}

	int internal_len = atoi(arr_size_str);

	if (matrix_str == 0)
		return 0;

	if (len != internal_len)
		Log.Warning("Requested matrix size doesn't match declaration in config file");

	int i = 0;
	for (; i < len; ++i) {
		SkipWhitespace(matrix_str);
		data[i] = f(matrix_str);
		SkipNonWhitespace(matrix_str);
	}
	return i;
}

int _Config::GetIntArray(char const * const name, int * pData, int const len) const {
	return ParseMatrix<int>(this, name, pData, len, atoi);
}

// strangly enough, the c standard declares atof as returning double
float _atof(char const * const str) {
	return (float) atof(str);
}

int _Config::GetFloatArray(char const * const name, float * pData, int const len) const {
	return ParseMatrix<float>(this, name, pData, len, _atof);
}

int _Config::GetDoubleArray(char const * const name, double * pData, int len) const {
	return ParseMatrix<double>(this, name, pData, len, atof);
}

bool _Config::Exists(char const * const name) const {
	return InternalExists(name);
}

bool _Config::HasArray(char const * const name) const {
	return config_arrays->find(std::string(name)) != config_arrays->end();
}

// Simple config file parser
void _Config::ParseFile(char const * const filename) {
	if (sealed) {
		Log.Error("Cant parse config file - config has been sealed.");
		return;
	}

	if (!FileExists(filename)) {
		Log.Error("Couldnt find config file \"%s\".", filename);
		return;
	}

	Log.Message("Parsing config file %s", filename);
	char const * config_text = 0;
	int map_id = CreateMapping(filename, config_text);

	while (*config_text != 0) {
		if (*config_text == '#')
			SkipLine(config_text, false);

		SkipWhitespace(config_text);

		if (*config_text != 0 && *config_text != '#') {
			char const * line = config_text;
			int n = 0;

			while (*(config_text + n) != 0 && *(config_text + n) != ' ')
				++n;

			std::string name = std::string(config_text, n);
			config_text += n;

			SkipWhitespace(config_text);

			if (*config_text == 0 || *config_text != '=') {
				SkipLine(config_text, false);
				int length = config_text - line - 1;
				Log.Warning("Ignoring malformed config entry in '%.*s'. Expected format: name=value.", length, line);
			} else {
				++config_text;

				SkipWhitespace(config_text);

				char const * start = config_text;

				SkipNonWhitespace(config_text);

				ulong len = (ulong) config_text - (ulong) start;
				std::string value(start, len);

				SkipLine(config_text, true);
				std::string arr_str;
				if (*config_text == '{') {
					++config_text;
					start = config_text;

					while (*config_text != '\0' && *config_text != '}')
						++config_text;

					if (*config_text != '}')
						Log.Error("'}' expected.");
					else {
						++config_text;
					}
					arr_str = std::string(start, config_text - start);
				}

				InternalAdd(name, value, arr_str);
				Log.Verbose("%s %s %s", name.c_str(), value.c_str(), arr_str.c_str());
			}
		}
	}

	CloseMapping(map_id);
}

// Initializes config by parsing arguments from command line
// and a config file
_Config::_Config(int argc, char * argv[]) {
	initialized = false;
	sealed = false;

	config_map = new TConfigMap();
	config_arrays = new TConfigMap();

	ParseArguments(argc, argv);

	if (InternalExists("config")) {
		ParseFile(InternalGet("config").c_str());
	}

	//Set default values
	Default("cpu_threads", 1);

	Default("kmer", 13);
	Default("kmer_skip", 2);
	Default("mode", 0);
//	Default("topn", 1);
	Default("paired", 0);
	Default("bs_mapping", 0);

	Default("max_insert_size", 1000);
	Default("min_insert_size", 0);
	Default("qry_start", 0);
	Default("qry_count", -1);

	Default("max_equal", 1);
	if (GetInt("bs_mapping") != 1) {
		Default("score_match", 5);
		Default("score_mismatch", -2);
		Default("score_gap_read", -6);
		Default("score_gap_ref", -6);
	} else {
		Log.Message("Using bs-mapping scoring scheme");
		Default("score_match", 4);
		Default("score_mismatch", -2);
		Default("score_gap_read", -10);
		Default("score_gap_ref", -10);
	}
	Default("score_match_tt", 4);
	Default("score_mismatch_tc", 4);

	//Silent
	Default("dualstrand", 1);
	Default("overwrite", 1);
	Default("ref_mode", -1);
//	Default("bowtie_mode", 0);

	//GPU
	Default("ocl_threads", 1);
	Default("block_multiplier", 2);
	Default("step_count", 4);

	//Filter
	Default("min_identity", 0);
	Default("min_residues", 0);

	//Output Input options
	Default("parse_all", 0);
	Default("hard_clip", 0);
	Default("silent_clip", 0);

	initialized = true;

	if (Exists("max_consecutive_indels")) {
		if (Exists("corridor")) {
			Log.Warning("Unable to use parameter 'max-consecutive-indels. Please remove the 'corridor' entry from the config file and use 'max-consecutive-indels'.");
		} else {
			std::stringstream ss;
			ss << GetInt("max_consecutive_indels", 5, 40) * 2;
			InternalAdd("corridor", ss.str(), "", true);
		}
	}
	if (Exists("bam")) {
		Default("format", 2);
	} else {
		Default("format", 1);
	}

	//Add full command line for CIGARWriter
	std::stringstream cmdLine;
	for (TConfigMap::iterator it = config_map->begin(); it != config_map->end(); ++it) {
		cmdLine << " --" << it->first + " " << it->second;
		if (config_arrays->count(it->first) != 0) {
			cmdLine << " {" << (*config_arrays)[it->first];
		}
	}
	InternalAdd("cmdline", cmdLine.str(), "", true);
	Log.Message("Parameter: %s", cmdLine.str().c_str());
}

//IConfig::~IConfig() {}

_Config::~_Config() {
	delete config_arrays;
	delete config_map;
}

void _Config::ParseArguments(int argc, char * argv[]) {

	if (argc == 1) {
		Help();
	}

	while (1) {
		int index = -1;
		struct option * opt = 0;
		int result = getopt_long(argc, argv, getopt_short, long_options, &index);
		if (result == '?') {
//			Log.Message("Unkown parameter %c", optopt);
			Help();
		} else if (result == ':') {
			Log.Message("Argument missing for parameter %c", optopt);
			Help();
		} else if (result != -1) {
			//Log.Message("switch %c, optind %d, opterr %d, optopt %d, arg %s", result, optind, opterr, optopt, optarg);
			if (result != 0) {
				//Short
				index = 0;
				while ((opt = (struct option *) &(long_options[index])) != 0 && opt->val != result) {
					index += 1;
				}
			}

			//Long
			opt = (struct option *) &(long_options[index]);
			if (opt->has_arg == required_argument) {
				if (optarg == 0 || *optarg == '-') {
					Log.Message("Argument missing for parameter %c", result);
					Help();
				}
				InternalAdd(opt->name, optarg, "", false);
			} else {
				//Hack for optional arguments without =
				if (strcmp(opt->name, "gpu") == 0) {
					optarg = argv[optind];

					if (optarg == 0 || *optarg == '-') {
						InternalAdd(opt->name, "1", " 0 }", false);
					} else {
						//Options that require an array
						int elementCount = 1;
						std::stringstream arrayData;
						for (uint i = 0; i < strlen(optarg); ++i) {
							if (optarg[i] == ',') {
								elementCount += 1;
								optarg[i] = ' ';
							}
						}
						arrayData << " " << optarg << " }";
						Log.Verbose("%s", arrayData.str().c_str());

						std::stringstream value;
						value << elementCount;
						InternalAdd(opt->name, value.str(), arrayData.str(), false);
						optind += 1;
					}
				} else {
					if (optarg == 0 || *optarg == '-') {
						InternalAdd(opt->name, "1", "", false);
					} else {
						InternalAdd(opt->name, optarg, "", false);
					}
				}
			}
		} else {
			break;
		}
	}

	/* print all other parameters */
	while (optind < argc) {
		printf("other parameter: <%s>\n", argv[optind++]);
	}

}

char const * const defaultConfigText =
		"################################################################################\n\
#                                                                              #\n\
# NGM Config Example                                                           #\n\
#                                                                              # \n\
# Lines starting with \"#\" are comments and thus ignored                        #\n\
#                                                                              #\n\
################################################################################\n\
# General                                                                      #\n\
################################################################################\n\
overwrite = 1\n\
gpu = 1 { 0 }\n\
mason_path = 	/software/ngm/ngm/mason/	\n\
\n\
\n";

