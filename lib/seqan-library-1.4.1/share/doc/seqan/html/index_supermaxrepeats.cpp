#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
	String<char> myString = "How many wood would a woodchuck chuck.";

	typedef Index< String<char> > TMyIndex;
	TMyIndex myIndex(myString);

	Iterator< TMyIndex, SuperMaxRepeats >::Type myRepeatIterator(myIndex, 3);

	while (!atEnd(myRepeatIterator)) 
	{
		for(unsigned i = 0; i < countOccurrences(myRepeatIterator); ++i)
			std::cout << getOccurrences(myRepeatIterator)[i] << ", ";

		std::cout << repLength(myRepeatIterator) << "   ";

		std::cout << "\t\"" << representative(myRepeatIterator) << '\"' << std::endl;

		++myRepeatIterator;
	}

	return 0;
}
