#include <iostream>
#include <seqan/find.h>

using namespace seqan;

template <typename TAlgorithm>
void printAllOccs(String<char>& haystack, 
				  String<char>& needle)
{
	Finder<String<char> > finder(haystack);
	Pattern<String<char>, TAlgorithm> pattern(needle);
	while (find(finder, pattern)) 
	{
		std::cout << position(finder) << ", ";
	}
	std::cout << std::endl;
}

int main() 
{
	String<char> haystack = "send more money!";
	String<char> needle = "mo";

	printAllOccs<Horspool>(haystack, needle);
	printAllOccs<BomAlgo> (haystack, needle);
	printAllOccs<BndmAlgo>(haystack, needle);
	printAllOccs<ShiftAnd>(haystack, needle);
	printAllOccs<ShiftOr> (haystack, needle);

	return 0;
}

