#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main() {
	String<unsigned int> seq;
	appendValue(seq, 5); appendValue(seq, 3); appendValue(seq, 4);
	appendValue(seq, 9); appendValue(seq, 6); appendValue(seq, 2);
	appendValue(seq, 1); appendValue(seq, 8); appendValue(seq, 7);
	appendValue(seq, 10);
	typedef Position<String<unsigned int> >::Type TPosition;
	String<TPosition, Block<> > pos;
	longestIncreasingSubsequence(seq,pos);
	for(int i = 0; i<(int) length(seq); ++i) {
		std::cout << seq[i] << ',';
	}
	std::cout << std::endl;
	std::cout << "Lis: " << std::endl;
	for(int i = length(pos)-1; i>=0; --i) {
		std::cout << seq[pos[i]] <<  ',';
	}
	std::cout << std::endl;
	return 0;
}
