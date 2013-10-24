#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
	StringSet< String<char> > mySet; 
	resize(mySet, 3); 
	mySet[0] = "tobeornottobe"; 
	mySet[1] = "thebeeonthecomb"; 
	mySet[2] = "beingjohnmalkovich"; 

	Index< StringSet<String<char> > > myIndex(mySet); 
	Finder< Index<StringSet<String<char> > > > myFinder(myIndex);

	std::cout << "hit at ";
	while (find(myFinder, "be")) 
		std::cout << position(myFinder) << "  ";
	std::cout << std::endl;
	
	return 0;
}
