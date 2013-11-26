#include <iostream>
#include <seqan/find.h>

using namespace seqan;

int main() 
{
	String<char> hayst = "If you must cross a course cross cow across a crowded cow crossing, "
						 "cross the cross coarse cow across the crowded cow crossing carefully.";
	String<char> ndl = "cr?o[uw]";
	Finder<String<char> > finder(hayst);
	Pattern<String<char>, WildShiftAnd> pattern(ndl);

	while (find(finder, pattern)) {
		std::cout << position(finder) << "\n";
	}
	return 0;
}

