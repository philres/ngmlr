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

#include "CS.h"
#include "Version.h"

#include "Timing.h"

#include "Debug.h"
#include "UpdateCheck.h"

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
		Log.Error("No output file (-o) specified. Writing to stdout.");
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


int main(int argc, char * argv[]) {

	std::stringstream version;
	version << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_BUILD;

	char const * arch[] = { "x86", "x64" };
	char const * build = (cDebug) ? " (DEBUG)" : "";
	Log.Message("NextGenMap-LR %s", version.str().c_str());
	Log.Message("Startup : %s%s (build %s %s)", arch[sizeof(void*) / 4 - 1], build, __DATE__, __TIME__);

	Log.Message("Starting time: %s", currentDateTime().c_str());

//try {
	Timer tmr;
	tmr.ST();

	_NGM::AppName = argv[0];

	InitPlatform();

// Initialization:
	try {
		_config = new _Config(argc, argv); // Parses command line & parameter file
		_log = &Log;

#ifdef DEBUGLOG
		if (Config.Exists(LOG)) {
			Log.Message("Writing debug log to stdout. Please use -o/--output for SAM/BAM output.");
			//Init checks if first parameter is != 0. Thus "LOG" is passed as a dummy string.
			_Log::Init("LOG", Config.GetInt(LOG_LVL));// Inits logging to file
		}
#else
		_Log::Init(0, 0); // Inits logging to file
#endif
	} catch (...) {
		Help();
	}

	if ( Config.GetInt("update_check")) {
		UpdateCheckInterface::remoteCheck();
		exit(0);
	}

	Log.setColor(Config.Exists("color"));

	if (!Config.Exists("qry") || CheckOutput()) {
		NGM; // Init Core

		NGM.InitProviders();

		if (!Config.Exists("qry") && !(Config.Exists("qry1") && Config.Exists("qry2"))) {
			Log.Message("Finished building hash table.");
			Log.Message("No qry file specified. If you want to map single-end data use -q/--qry. If you want to map paired-end data, either use -q/--qry and -p or --qry1 and --qry2");
		} else {
			try {
				NGM.StartThreads();

				bool const progress = Config.GetInt(PROGRESS) == 1;

				if(!progress) {
					Log.Message("No progress information (use --progress)");
				}
				Log.Message("Mapping reads...");

				NGM.MainLoop();

//				if(NGM.Stats->invalidAligmentCount > 0) {
//					Log.Message("Unsuccessful alignment computation %d", NGM.Stats->invalidAligmentCount);
//				}

				int discarded = NGM.GetReadReadCount() - (NGM.GetMappedReadCount()+NGM.GetUnmappedReadCount());
				if (discarded != 0) {
					Log.Message("Done (%i reads mapped (%.2f%%), %i reads not mapped (%i discarded), %i lines written)(elapsed: %ds, %d r/s)", NGM.GetMappedReadCount(), (float)NGM.GetMappedReadCount() * 100.0f / (float)(std::max(1, NGM.GetMappedReadCount()+NGM.GetUnmappedReadCount() + discarded)),NGM.GetUnmappedReadCount() + discarded, discarded, NGM.GetWrittenReadCount(), (int) tmr.ET(), (int)(NGM.GetMappedReadCount() * 1.0f / tmr.ET()));
				} else {
					Log.Message("Done (%i reads mapped (%.2f%%), %i reads not mapped, %i lines written)(elapsed: %ds, %d r/s)", NGM.GetMappedReadCount(), (float)NGM.GetMappedReadCount() * 100.0f / (float)(std::max(1, NGM.GetMappedReadCount()+NGM.GetUnmappedReadCount())),NGM.GetUnmappedReadCount(),NGM.GetWrittenReadCount(), (int) tmr.ET(), (int)(NGM.GetMappedReadCount() * 1.0f / tmr.ET()));
				}

				NGM.StopThreads();

			} catch (...) {
				Log.Error("Unhandled exception in control thread.");
			}
		}
	}

	if (! Config.GetInt("update_check"))
		UpdateCheckInterface::reminder();

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
  ngm -q <reads> -r <reference> -o <output> [parameter]\n\
\n\
Input:\n\
\n\
 -r/--reference <path>         Path to the reference genome\n\
                               (format: FASTA, can be gzipped).\n\
 -q/--qry  <path>              Path to the read file. If the file contains\n\
                               interleaved mates use -p/--paired.\n\
                               find the max. read length. (default: off)\n\
\n\
Output:\n\
\n\
 -o/--output <path>            Path to output file.\n\
 -b/--bam                      Output BAM instead of SAM.\n\
 --keep-tags                   Copy BAM/SAM tags present in input file to\n\
                               output file (default: off)\n\
 --no-unal                     Don't print unaligned reads to output file.\n\
 --hard-clip                   Hard instead of soft clipping in SAM output\n\
 --silent-clip                 Hard clip reads but don't add clipping\n\
                               information to CIGAR string\n\
 --rg-id <string>              Adds RG:Z:<string> to all alignments in SAM/BAM \n\
 --rg-sm <string>              RG header: Sample\n\
 --rg-lb <string>              RG header: Library\n\
 --rg-pl <string>              RG header: Platform\n\
 --rg-ds <string>              RG header: Description\n\
 --rg-dt <string>              RG header: Date (format: YYYY-MM-DD)\n\
 --rg-pu <string>              RG header: Platform unit\n\
 --rg-pi <string>              RG header: Median insert size\n\
 --rg-pg <string>              RG header: Programs\n\
 --rg-cn <string>              RG header: sequencing center\n\
 --rg-fo <string>              RG header: Flow order\n\
 --rg-ks <string>              RG header: Key sequence\n\
\n\
General:\n\n\
 -t/--threads <int>            Number of candidate search threads\n\
 -i/--min-identity <0-1>       All reads mapped with an identity lower than\n\
                               this threshold will be reported as unmapped\n\
                               (default: 0.65)\n\
 -R/--min-residues <int/float> All reads mapped with lower than\n\
                               <int> or <float> * read length residues\n\
                               will be reported as unmapped. (default: 0.5)\n\
 -Q/--min-mq <int>             All reads mapped with lower than <int>\n\
                               mapping quality will be reported as unmapped.\n\
                               (default: 0)\n\
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
 --bin-size <n>                Sets the size of the grid NextgenMap uses\n\
                               during CMR search to: 2^n (default: 2) \n\
\n\
Other:\n\
\n\
 --update-check                Perform an online check for a newer version of NGM\n\
 --color                       Colored text output (default: off)\n\
 --progress                    Print progress info while mapping\n\
                               (default: off)\n\
\n\
\n";
	fprintf(stderr, "%s", help_msg);
	exit(1);
}

// actually platform specific.../care
uloc const FileSize(char const * const filename) {
	FILE * fp = fopen(filename, "rb");
	if (fp == 0) {
		Log.Warning("Tried to get size of nonexistent file %s", filename);
		return 0;
	}

	if (fseek(fp, 0, SEEK_END) != 0)
		return 0;

#ifdef __APPLE__
	uloc end = ftello(fp);
#else
	uloc end = ftello64(fp);
#endif

	fclose(fp);
	return end;
}
