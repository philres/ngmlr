#include <stdlib.h>
#include <stdio.h>

#ifdef _WIN32
#include <conio.h>
#endif

#ifdef _WIN32
#include <crtdbg.h>
#endif

#include <unistd.h>
#include <sstream>

#include "NGM.h"

#include "SW.h"
#include "CS.h"
#include "Output.h"
#include "Version.h"

#include "Timing.h"

#include "Debug.h"

#undef module_name
#define module_name "MAIN"

bool CheckOutput() {
	if (Config.Exists("output")) {
		if (!(Config.Exists("overwrite") && (Config.GetInt("overwrite") == 1))) {
			if (FileExists(Config.GetString("output"))) {
				Log.Error("Output file already exists and overwrite set to false.");
				return false;
			}
		}
	} else {
		Log.Error("No output file (-o) specified.");
		Help();
		return false;
	}
	if (Config.Exists("qry")) {
		if (!FileExists(Config.GetString("qry"))) {
			Log.Error("Query file (%s) does not exist.", Config.GetString("qry"));
			Help();
			return false;
		}
	}
	return true;
}

#ifdef NDEBUG
bool cDebug = false;
#else
bool cDebug = true;
#endif

ILog const * _log = 0;
IConfig * _config = 0;

//char const * const version = "0.4.4";

int main(int argc, char * argv[]) {

	std::stringstream version;
	version << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_BUILD;

	char const * arch[] = { "x86", "x64" };
	char const * build = (cDebug) ? " (DEBUG)" : "";
	Log.Message("NextGenMap %s", version.str().c_str());
	Log.Message("Startup : %s%s (build %s %s)", arch[sizeof(void*) / 4 - 1], build, __DATE__, __TIME__);

	//try {
	Timer tmr;
	tmr.ST();

	_NGM::AppName = argv[0];

	InitPlatform();

	// Initialization:
	_config = new _Config(argc, argv); // Parses command line & parameter file
	_log = &Log;
	_Log::Init(); // Inits logging to file

//	//Restart NGM and set proper environment variables
//	if (Config.Exists("mason_path") && !Config.Exists("skip_env")) {
//		char * newArgv[argc + 3];
//		for (int i = 0; i < argc; ++i) {
//			newArgv[i] = argv[i];
//		}
//		newArgv[argc] = "--skip-env";
//		newArgv[argc + 1] = "1";
//		newArgv[argc + 2] = (char *) 0;
//
//		int const length = 3;
//		char * envParms[length] = { NULL, NULL };
//
//		envParms[0] = new char[10000];
//		envParms[1] = new char[10000];
//		int index = 0;
//
//		index += sprintf(envParms[0], "LD_LIBRARY_PATH=");
//		if (Config.Exists("mason_path")) {
//			index += sprintf(envParms[0] + index, "%s", Config.GetString("mason_path"));
//		}
//		if (Config.Exists("lib_path")) {
//			index += sprintf(envParms[0] + index, ":%s", Config.GetString("lib_path"));
//		}
//		if (Config.Exists("vendor_path")) {
//			sprintf(envParms[1], "OPENCL_VENDOR_PATH=%s", Config.GetString("vendor_path"));
//			envParms[2] = NULL;
//		} else {
//			envParms[1] = NULL;
//		}
//
////			for (int i = 0; i < length; ++i) {
////				Log.Message("ENV: %s", envParms[i]);
////			}
//		//BrInitError error;
//		//br_init(&error);
//		//Log.Message("Restarting with environment variables (%s).", br_find_exe(argv[0]));
//		Log.Message("Restarting with environment variables (%s).", argv[0]);
//		int ret = execve(newArgv[0], newArgv, envParms);
//		Log.Error("Couldn't restart NGM. Make sure that all environment variables are set correctly (%d).", ret);
//	}

	if (Config.Exists("master_cpu"))
		NGMSetThreadAffinity(0, Config.GetInt("master_cpu"));

	if (!Config.Exists("qry") || CheckOutput()) {
		NGM; // Init Core

		NGM.InitProviders();

		if (!Config.Exists("qry")) {
			Log.Message("Finished building hash table. No qry file specified. Exiting.");
		} else {
			Log.Message("Core initialization complete");
			NGM.StartThreads();

			NGM.MainLoop();
			Log.Message("Scores computed: %ld", SW::scoreCount);
			Log.Message("Alignments computed: %ld", Output::alignmentCount);
#ifdef INSTANCE_COUNTING
			Log.Green("Counts:");
			Log.Message("MappedRead count = %i (max %i)", MappedRead::sInstanceCount, MappedRead::maxSeqCount);
			Log.Message("LocationScore count = %i", LocationScore::sInstanceCount);
#endif
			int discarded = NGM.GetReadReadCount() - NGM.GetWrittenReadCount();
			if (discarded != 0) {
				Log.Warning("Reads discarded: %d", discarded);
			}

			Log.Message("Done (%i reads mapped (%.2f%%), %i reads not mapped, %i reads written)(elapsed: %fs)", NGM.GetMappedReadCount(), (float)NGM.GetMappedReadCount() * 100.0f / (float)(std::max(1, NGM.GetMappedReadCount()+NGM.GetUnmappedReadCount())),NGM.GetUnmappedReadCount(), NGM.GetWrittenReadCount(), tmr.ET());
			//} catch (...) {
			//	Log.Error("Unhandled exception in control thread");
			//}
		}
	}
	CS::Cleanup();
	_SequenceProvider::Cleanup();
	delete _config;
	delete _NGM::pInstance;
	_Log::Cleanup();

	CleanupPlatform();

	return 0;
}

void Help() {
	char const * const help_msg =
			"\
\
Usage:\
  ngm [-c <path>] -q <reads> -r <reference> -o <output> [parameter]\n\
\n\
Input/Output:\n\
\n\
  -c/--config <path>            Path to the config file. The config file contains all advanced options.\n\
                                If this parameter is omitted, default values will be used.\n\
  -r/--reference <path>         Path to the reference genome (format: FASTA, can be gzipped).\n\
  -q/--query <path>             Path to the read file. Valid input formats are: FASTA/Q (gzipped), SAM/BAM\n\
                                If the query file is omitted, NGM will only preprocess the reference.\n\
  -p/--paired                   Input data is paired end. (default: off)\n\
  -I/--min-insert-size          The minimum insert size for paired end alignments (default: 0)\n\
  -X/--max-insert-size          The maximum insert size for paired end alignments (default: 1000)\n\
\n\
Output:\n\
\n\
  -o/--output <path>            Path to output file.\n\
  -b/--bam						Output BAM instead of SAM.\n\
  --hard-clip                   Use hard instead of soft clipping for SAM output\n\
  --silent-clip                 Hard clip reads but don't add clipping information to CIGAR string\n\
\n\
General:\n\n\
  -t/--threads <int>            Number of candidate search threads\n\
  -s/--sensitivity <float>      A value between 0 and 1 that determines the number of possible mapping\n\
                                locations that will be evaluated with an sequence alignment.\n\
                                0 means that all possible mapping locations will be evaluated\n\
                                1 means that only the best possible mapping location(s) will be evaluated\n\
                                Higher values will reduce the runtime but also have a negativ effect on mapping sensitivity.\n\
                                (default: estimated from input data)\n\
  -i/--min-identity <0-1>       All reads mapped with an identity lower than this threshold will be reported as unmapped (default: 0.75)\n\
  -R/--min-residues <int>       All reads mapped with lower than <int> residues will be reported as unmapped (default: 0.0)\n\
  -g/--gpu [int,...]            Use GPU(s) for alignment computation (GPU Ids are zero-bassed!).\n\
                                   With -g or --gpu GPU 0 will be used.\n\
                                   With -g 1 or --gpu 1 GPU 1 will be used.\n\
                                   With -g 0,1 or --gpu 0,1 GPU 0 and 1 will be used.\n\
                                If -g/--gpu is ommitted, alignments will be computed on the CPU\n\
  --bs-mapping                  Enables bisulfite mapping (For bs-mapping kmer-skip will be applied to\n\
                                the reads instead of the reference sequence).\n\
  --bs-cutoff <int>             Max. number of Ts in a k-mer. All k-mers were the number of Ts is higher than <int> are ignored (default: 8)\n\
  -h/--help                     Prints help and aborts program\n\n\
Advanced settings:\n\
\n\
  -k/--kmer [10-14]             Kmer length in bases. (default: 13)\n\
  --kmer-skip <int>             Number of k-mers to skip when building the lookup table from\n\
                                the reference(default: 2)\n\
  -m/mode [0|1]                 Alignment mode: 0 = local, 1 = semi-global. (default: 0)\n\
  -C/--max-consec-indels <int>  Maximum number of consecutive indels allowed. (default: computed from input)\n\
\n\
\n";
	printf(help_msg);
	exit(0);
}

// actually platform specific.../care
ulong const FileSize(char const * const filename) {
	FILE * fp = fopen(filename, "rb");
	if (fp == 0) {
		Log.Warning("Tried to get size of nonexistant file %s", filename);
		return 0;
	}

	if (fseek(fp, 0, SEEK_END) != 0)
		return 0;

	ulong end = ftell(fp);
	fclose(fp);
	return end;
}
