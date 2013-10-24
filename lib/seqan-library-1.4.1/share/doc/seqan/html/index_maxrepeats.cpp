#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
	String<char> myString = "How many wood would a woodchuck chuck.";

	typedef Index< String<char> > TMyIndex;
	TMyIndex myIndex(myString);

	typedef Iterator< TMyIndex, MaxRepeats >::Type TMaxRepeatIterator;
	TMaxRepeatIterator myRepeatIterator(myIndex, 3);

	while (!atEnd(myRepeatIterator)) 
	{
		Iterator<TMaxRepeatIterator>::Type myRepeatPair(myRepeatIterator);
		while (!atEnd(myRepeatPair)) {
			std::cout << *myRepeatPair << ", ";
			++myRepeatPair;
		}

		std::cout << repLength(myRepeatIterator) << "   ";

		std::cout << "\t\"" << representative(myRepeatIterator) << '\"' << std::endl;

		++myRepeatIterator;
	}

	return 0;
}
