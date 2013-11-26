#include <iostream>
#include <fstream>
#include <cstdio>

#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
	std::FILE * fl = std::fopen("testfile.fa", "wb");
    write(fl, "aacagtattagaccactaggaccct", "a test file", Fasta());
	close (fl);
	std::fstream fstrm;
	fstrm.open("testfile.fa", std::ios_base::in | std::ios_base::binary);
	String<char> fasta_tag;
	String<Dna> fasta_seq;
	readMeta(fstrm, fasta_tag, Fasta());
	std::cout << fasta_tag << "\n";	//prints "a test file"
	read(fstrm, fasta_seq, Fasta());
	std::cout << fasta_seq << "\n";	//prints the sequence
	fstrm.close();
	String<Dna, FileReader<Fasta> > fr("testfile.fa");
	std::cout << fr << "\n";			//prints the sequence
	return 0;
}

