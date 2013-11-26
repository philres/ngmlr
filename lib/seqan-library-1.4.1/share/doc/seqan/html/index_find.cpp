#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
	Index< String<char> > index_esa("tobeornottobe");
	Finder< Index< String<char> > > finder_esa(index_esa);

	std::cout << "hit at ";
	while (find(finder_esa, "be"))
		std::cout << position(finder_esa) << " ";
	std::cout << std::endl;

	typedef Index< String<char>, IndexQGram< UngappedShape<2> > > TQGramIndex;
	TQGramIndex index_2gram("tobeornottobe");
	Finder< TQGramIndex > finder_2gram(index_2gram);

	std::cout << "hit at ";
	while (find(finder_2gram, "be"))
		std::cout << position(finder_2gram) << " ";
	std::cout << std::endl;

	return 0;
}
