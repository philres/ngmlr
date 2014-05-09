/*
 * main.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: philipp_
 */
#include "interleave-pairs.h"

#include <zlib.h>
#include <stdio.h>
#include <limits.h>
#include <iostream>
#include <map>
#include <cstring>
#include <tclap/CmdLine.h>

#include "Log.h"
#include "IParser.h"
#include "PrefixTable.h"
#include "SequenceProvider.h"
#include "CS.h"


using std::cout;
using std::endl;
using std::map;
using std::string;

#undef module_name
#define module_name "MAIN"

struct nString {
	char * s;
	size_t n;
};

struct Sequence {
	nString name, seq, qlty;
	int fragPos;
};

IRefProvider * m_RefProvider = 0;

FILE * ofp;

void CountKmer(ulong prefix, uint pos, ulong mutateFrom, ulong mutateTo, void* data) {



	RefEntry * m_entry = new RefEntry(0);
	m_entry->nextEntry = new RefEntry(0);

	RefEntry const * cur = m_RefProvider->GetRefEntry(prefix, m_entry);

	SequenceLocation loc;
	loc.m_Location = pos;
	if(!SequenceProvider.convert(loc)) {
		Log.Message("Couldn't decode sequence location.");
		Fatal();
	}

	int refNameLength = 0;
	char const * refName = SequenceProvider.GetRefName(loc.getrefId(), refNameLength);
	fprintf(ofp, "%s\t%u\t%u\t%d\t%d\n", refName, loc.m_Location, loc.m_Location + 1, cur->refCount, cur->nextEntry->refCount);
	//	int * freq = (int *) data;
//	if (prefix == lastPrefix) {
//		int currentBin = GetBin(pos);
//		if (currentBin != lastBin || lastBin == -1) {
//			freq[prefix] += 1;
//		} else {
//			skipCount += 1;
////			Log.Message("Prefix %d (skip):\t%d (%d)\t%d (%d)\t(%ld)", prefix, lastPos, lastBin, pos, currentBin, skipCount);
//		}
//		lastBin = currentBin;
//		lastPos = pos;
//	} else {
//		lastBin = -1;
//		lastPos = -1;
//		freq[prefix] += 1;
//	}
//	lastPrefix = prefix;
}


int kmer_distribution(int argc, char **argv) {

	try {

		TCLAP::CmdLine cmd("Interleaves paired end reads from two FASTA/Q files into one FASTQ file.", ' ', "0.1", false);

		//TCLAP::ValueArg<std::string> leftArg("1", "m1", "Upstream mates (FASTA/Q)", true, "", "file");
		//TCLAP::ValueArg<std::string> rightArg("2", "m2", "Downstream mates (FASTA/Q)", true, "", "file");

		TCLAP::ValueArg<std::string> refArg("r", "ref", "Reference file", true, "", "file");
		//TCLAP::ValueArg<std::string> unpairedArg("u", "unpaired", "Write reads without mate to this file.", false, "", "file");

		TCLAP::ValueArg<std::string> outArg("o", "out", "Output file", true, "", "file");

		TCLAP::ValueArg<int> kmerArg("k", "kmer", "K-mer length", false, 13,
				"int");

		//TCLAP::SwitchArg noprogressArg("", "noprogress", "Suppress progress output.", cmd, false);

		//TCLAP::SwitchArg forceArg("f", "force", "Force finishing even if no pairs are found.", cmd, false);

		//cmd.add(delimiterArg);
		//cmd.add(unpairedArg);
		//cmd.add(outArg);
		//cmd.add(rightArg);
		cmd.add(refArg);
		cmd.add(outArg);
		cmd.add(kmerArg);

		cmd.parse(argc, argv);

		_log = &Log;
		_Log::Init(0,0); // Inits logging to file

		((_Config*) _config)->Default("bs_mapping", 0);
		((_Config*) _config)->Default("skip_save", 0);
		((_Config*) _config)->Default("kmer_maxFreq", 1000000);

		//((_Config*) _config)->Default("ref", "/project/ngs-work/meta/reference/genomes/eck12_MG1655_ecoli/eck12_MG1655_ecoli.fasta");
		((_Config*) _config)->Default("ref", refArg.getValue().c_str());


		((_Config*) _config)->Override("kmer", kmerArg.getValue());
		((_Config*) _config)->Override("kmer_skip", 0);


		CS::prefixBasecount = Config.GetInt("kmer", 4, 32);
		CS::prefixBits = CS::prefixBasecount * 2;
		CS::prefixMask = ((ulong) 1 << CS::prefixBits) - 1;

		SequenceProvider.Init(); // Prepares input data

//		char const * refName = 0;
//		int refNameLength = 0;
//
//		for (int i = 0; i < SequenceProvider.GetRefCount(); ++i) {
//			refName = SequenceProvider.GetRefName(i, refNameLength);
//			Log.Message("@SQ\tSN:%.*s\tLN:%d\n", refNameLength, refName, SequenceProvider.GetRefLen(i));
//			++i;
//		}
//
//		char * test = new char[100];
//		SequenceProvider.DecodeRefSequence(test, 0, 10000, 100);
//		Log.Message("%s", test);

		m_RefProvider = new CompactPrefixTable(true, false);



		int length = 100;
		int * freq = new int[length];
		memset(freq, 0, length);


		Log.Message("Processing: ");
		ofp = fopen(outArg.getValue().c_str(), "w");
		for (int i = 0; i < SequenceProvider.GetRefCount(); i = i + 2) {
//			lastPrefix = 111111;
//			lastBin = -1;
//			m_CurGenSeq = i;
//
//			if (!DualStrand || !(m_CurGenSeq % 2)) {
//
				uint offset = SequenceProvider.GetRefStart(i);
				uint len = SequenceProvider.GetRefLen(i);
				int m_RefSkip = 0;
				char * seq = new char[len + 2];
				SequenceProvider.DecodeRefSequence(seq, i, offset, len);

				int refNameLength = 0;
				char const * refName = SequenceProvider.GetRefName(i, refNameLength);
				Log.Message("%s", refName);
//
				CS::PrefixIteration(seq, len, &CountKmer, 0, 0, freq, m_RefSkip, offset);
//				delete[] seq;
//				seq = 0;
//			}
		}
//		return freq;
		delete[] freq;

		fclose(ofp);
//

	} catch (TCLAP::ArgException &e) {
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	} catch (std::ios_base::failure &e) {
		std::cerr << "Error: " << e.what() << std::endl;
	} catch(char const * msg) {
		std::cerr << "Exception: " << msg << std::endl;
	}

	return 0;
}
