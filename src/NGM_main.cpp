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

//#include "SW.h"
#include "CS.h"
//#include "Output.h"
#include "Version.h"

#include "Timing.h"

#include "Debug.h"

#undef module_name
#define module_name "MAIN"

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	// Visit http://www.cplusplus.com/reference/clibrary/ctime/strftime/
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

	return buf;
}

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
	if((Config.Exists("qry1") && Config.Exists("qry2"))) {
		if (!FileExists(Config.GetString("qry1"))) {
			Log.Error("Query file (%s) does not exist.", Config.GetString("qry1"));
			Help();
			return false;
		}
		if (!FileExists(Config.GetString("qry2"))) {
			Log.Error("Query file (%s) does not exist.", Config.GetString("qry2"));
			Help();
			return false;
		}
	} else {
		if (Config.Exists("qry")) {
			if (!FileExists(Config.GetString("qry"))) {
				Log.Error("Query file (%s) does not exist.", Config.GetString("qry"));
				Help();
				return false;
			}
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

	Log.Message("Starting time: %s", currentDateTime().c_str());


	//try {
	Timer tmr;
	tmr.ST();

	_NGM::AppName = argv[0];

	InitPlatform();

	// Initialization:
	_config = new _Config(argc, argv); // Parses command line & parameter file
	_log = &Log;
	_Log::Init(); // Inits logging to file

	Log.setColor(Config.Exists("color"));

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

		if (!Config.Exists("qry") && !(Config.Exists("qry1") && Config.Exists("qry2"))) {
			Log.Message("Finished building hash table.");
			Log.Message("No qry file specified. If you want to map single-end data use -q/--qry. If you want to map paired-end data, either use -q/--qry and -p or --qry1 and --qry2.");
		} else {
			try {
				Log.Message("Core initialization complete");
				NGM.StartThreads();

				NGM.MainLoop();
				//Log.Message("Scores computed: %ld", SW::scoreCount);
				Log.Message("Alignments computed: %ld", AlignmentBuffer::alignmentCount);
#ifdef INSTANCE_COUNTING
				Log.Green("Counts:");
				Log.Message("MappedRead count = %i (max %i)", MappedRead::sInstanceCount, MappedRead::maxSeqCount);
				Log.Message("LocationScore count = %i", LocationScore::sInstanceCount);
#endif
				int discarded = NGM.GetReadReadCount() - (NGM.GetMappedReadCount()+NGM.GetUnmappedReadCount());
				if (discarded != 0) {
					Log.Warning("Reads discarded: %d", discarded);
				}

				Log.Message("Done (%i reads mapped (%.2f%%), %i reads not mapped, %i lines written)(elapsed: %fs)", NGM.GetMappedReadCount(), (float)NGM.GetMappedReadCount() * 100.0f / (float)(std::max(1, NGM.GetMappedReadCount()+NGM.GetUnmappedReadCount())),NGM.GetUnmappedReadCount(), NGM.GetWrittenReadCount(), tmr.ET());
			} catch (...) {
				Log.Error("Unhandled exception in control thread");
			}
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
  ngm [-c <path>] {-q <reads> [-p] | -1 <mates 1> -2 <mates 2>} -r <reference> -o <output> [parameter]\n\
\n\
Input/Output:\n\
\n\
 -c/--config <path>            Path to the config file. The config file contains\n\
                               all advanced options. If this parameter is\n\
                               omitted, default values will be used.\n\
 -r/--reference <path>         Path to the reference genome\n\
                               (format: FASTA, can be gzipped).\n\
 -q/--qry  <path>              Path to the read file. If the file contains\n\
                               interleaved mates use -p/--paired.\n\
 -1/--qry1 <path>              Path to the read file containing mates 1.\n\
 -2/--qry2 <path>              Path to the read file containing mates 2.\n\
                               Valid input formats are: FASTA/Q (gzipped),\n\
                               SAM/BAM. If the query file(s) is/are omitted,\n\
                               NGM will only pre-process the reference.\n\
 -p/--paired                   Input data is paired end.\n\
                               NOT required if -1/-2 are used. (default: off)\n\
 -I/--min-insert-size          The min insert size for paired end alignments\n\
                               (default: 0)\n\
 -X/--max-insert-size          The max insert size for paired end alignments\n\
                               (default: 1000)\n\
\n\
Output:\n\
\n\
 -o/--output <path>            Path to output file.\n\
 -b/--bam                      Output BAM instead of SAM.\n\
 --no-unal                     Don't print unaligned reads to output file.\n\
 --hard-clip                   Hard instead of soft clipping in SAM output\n\
 --silent-clip                 Hard clip reads but don't add clipping\n\
                               information to CIGAR string\n\
 -n/--topn                     Prints the <n> best alignments sorted by\n\
                               alignment score (default: 1)\n\
 --strata                      Only  output  the  highest  scoring  mappings\n\
                               for any  given  read,  up  to <n> mappings per\n\
                               read. If a read has more than <n> mappings with\n\
                               the same score, it is discarded and reported\n\
                               as unmapped.\n\
\n\
General:\n\n\
 -t/--threads <int>            Number of candidate search threads\n\
 -s/--sensitivity <float>      A value between 0 and 1 that determines the\n\
                               number of candidate mapping regions that will\n\
                               be evaluated with an sequence alignment.\n\
                                 0: all CMR(s) will be evaluated\n\
                                 1: only the best CMR(s) will be evaluated\n\
                               Higher values will reduce the runtime but also\n\
                               have a negative effect on mapping sensitivity.\n\
                               (default: estimated from input data)\n\
 -i/--min-identity <0-1>       All reads mapped with an identity lower than\n\
                               this threshold will be reported as unmapped\n\
                               (default: 0.65)\n\
 -R/--min-residues <int/float> All reads mapped with lower than\n\
                               <int> or <float> * read length residues\n\
                               will be reported as unmapped. (default: 0.5)\n\
 -g/--gpu [int,...]            Use GPU(s) for alignment computation\n\
                               NOTE: GPU Ids are zero-based!\n\
                                  -g or --gpu to use GPU\n\
                                  -g 1 or --gpu 1 to use GPU 1\n\
                                  -g 0,1 or --gpu 0,1 to use GPU 0 and 1\n\
                               If -g/--gpu is omitted, alignments will be\n\
                               computed on the CPU (default)\n\
 --bs-mapping                  Enables bisulfite mapping.\n\
                               For bs-mapping, kmer-skip will be applied to\n\
                               the reads instead of the reference sequence.\n\
 --bs-cutoff <int>             Max. number of Ts in a k-mer. All k-mers were\n\
                               the number of Ts is higher than <int> are\n\
                               ignored (default: 8)\n\
 -h/--help                     Prints help and aborts program\n\n\
\n\
Advanced settings:\n\
\n\
 -k/--kmer [10-14]             Kmer length in bases. (default: 13)\n\
 --kmer-skip <int>             Number of k-mers to skip when building the\n\
                               lookup table from the reference(default: 2)\n\
 --kmer-min <int>              Minimal number of k-mer hits required to\n\
                               consider a region a CMR. (default: 0)\n\
 --max-cmrs <int>              Reads that have more than <int> CMRs are\n\
                               ignored. (default: infinite)\n\
 --skip-save                   Don't save index to disk. Saves disk space but\n\
                               increases runtime for larger genomes.\n\
 -l/--local                    Perform local alignment. Ends might get clipped.\n\
                               (default: on)\n\
 -e/--end-to-end               Perform local alignment. Ends might get clipped.\n\
			                   (default: off)\n\
 --affine                      Use alignment algorithms that support affine gap\n\
                               costs. This will give you better alignments for\n\
                               longer indels but will also increase the runtime.\n\
                               (default: off)\n\
 -C/--max-consec-indels <int>  Maximum number of consecutive indels allowed.\n\
                               (default: computed based on avg. read length)\n\
 --fast-pairing                Mates are mapped individually. If the best\n\
                               alignments for the mates give a proper pair\n\
                               they are marked as paired in the output.\n\
                               If not they are reported as broken pair.\n\
 --pair-score-cutoff <0-1>     To find the best pairing all alignments of a\n\
                               read that have a score in the range of: \n\
                                 [top score * pair-score-cutoff; top score]\n\
                               are taken into account.\n\
 --match-bonus <int>           Match score\n\
                               (default: 10, affine: 10, bs-mapping: 4)\n\
 --mismatch penalty <int>      Mismatch penalty\n\
                               (default: 15, affine: 15, bs-mapping: 2)\n\
 --gap-read-penalty <int>      Penalty for a gap in the read\n\
                               (default: 20, affine: 33, bs-mapping: 10)\n\
 --gap-ref-penalty <int>       Penalty for a gap in the reference\n\
                               (default: 20, affine: 33, bs-mapping: 10)\n\
 --gap-extend-penalty <int>    Penalty for extending a gap\n\
			                   (default: 20, affine: 3, bs-mapping: 10)\n\
 --match-bonus-tt <int>        Only for bs-mapping (default: 4)\n\
 --match-bonus-tc <int>        Only for bs-mapping (default: 4)\n\
 --color                       Colored text output (default: off)\n\
 --no-progress                 Don't print progress info while mapping\n\
                               (default: off)\n\
\n\
\n";
	printf(help_msg);
	exit(0);
}

// actually platform specific.../care
ulong const FileSize(char const * const filename) {
	FILE * fp = fopen(filename, "rb");
	if (fp == 0) {
		Log.Warning("Tried to get size of nonexistent file %s", filename);
		return 0;
	}

	if (fseek(fp, 0, SEEK_END) != 0)
		return 0;

	ulong end = ftell(fp);
	fclose(fp);
	return end;
}
