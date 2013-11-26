#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
	StringSet< String<char> > mySet;
	resize(mySet, 3);
	mySet[0] = "SeqAn is a library for sequence analysis.";
	mySet[1] = "The String class is the fundamental sequence type in SeqAn.";
	mySet[2] = "Subsequences can be handled with SeqAn's Segment class.";

	typedef Index< StringSet<String<char> > > TMyIndex;
	TMyIndex myIndex(mySet);

	Iterator< TMyIndex, Mums >::Type myMUMiterator(myIndex, 3);
	String< SAValue<TMyIndex>::Type > occs;

	while (!atEnd(myMUMiterator)) {
		occs = getOccurrences(myMUMiterator);
		orderOccurrences(occs);
		
		for(unsigned i = 0; i < length(occs); ++i)
			std::cout << getValueI2(occs[i]) << ", ";

		std::cout << repLength(myMUMiterator) << "   ";

		std::cout << "\t\"" << representative(myMUMiterator) << '\"' << std::endl;

		++myMUMiterator;
	}

	return 0;
}
