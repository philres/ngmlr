#include <iostream>
#include <seqan/find.h>

using namespace seqan;

int main() 
{
	String<char> haystk("AACTTAACCTAA");
	String<char> ndl("CCT");

	Finder<String<char> > fnd(haystk);
	Pattern<String<char>, MyersUkkonen> pat(ndl);
	setScoreLimit(pat, -1);
	while (find(fnd, pat)) {
		std::cout << position(fnd) << ": " << getScore(pat) << "\n";
	}

	String<char> t = "babybanana";
	String<char> p = "babana";
	Finder<String<char> > finder(t);
	Pattern<String<char>, Myers<FindInfix> > pattern(p);
	while (find(finder, pattern, -2)) {
		std::cout << "end: " << endPosition(finder) << std::endl;
		while (findBegin(finder, pattern, getScore(pattern))) {
			std::cout << "begin: " << beginPosition(finder) << std::endl;
			std::cout << infix(finder) << " matches with score ";
			std::cout << getBeginScore(pattern) << std::endl;
		}
	}
	return 0;
}

