#include <iostream>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int main()
{
	typedef String<char> TText;
	TText str = "abcdefg";
	Iterator<TText, Rooted>::Type it = begin(str);
	it = begin(str, Rooted());
	std::cout << container(it);          //output: "abcdefg"
	goNext(it);
	std::cout << position(it);           //output: 7
	return 0;
}
