#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

int main()
{
	seqan::String<char> str = "admn";
	seqan::Iterator<seqan::String<char> >::Type it = begin(str);
	seqan::Iterator<seqan::String<char> >::Type itEnd = end(str);
	while (it != itEnd) {
		std::cout << *it;
		++it;
	}
	std::cout << std::endl;
	seqan::Iterator<seqan::String<char>, seqan::Rooted >::Type it2 = begin(str);
	for (goBegin(it2); !atEnd(it2); goNext(it2)) 
	{
		++value(it2);
	}
	goEnd(it2);
	while (!atBegin(it2))              
	{
		goPrevious(it2);
		std::cout << getValue(it2);
	}
	std::cout << std::endl;
	assignValue(begin(str), 'X');
	std::cout << str << std::endl;
	
	return 0;
}
