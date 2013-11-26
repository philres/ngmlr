#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
	String<char> text = "hello world!";
	String<char> pattern = "l";
	String<unsigned> sa;

	resize(sa, length(text));
	createSuffixArray(sa, text, Skew7());

	Pair<unsigned> hitRange;
	hitRange = equalRangeSA(text, sa, pattern);

	for(unsigned i = hitRange.i1; i < hitRange.i2; ++i)
		std::cout << sa[i] << " ";
	std::cout << std::endl;
 
	return 0;
}
