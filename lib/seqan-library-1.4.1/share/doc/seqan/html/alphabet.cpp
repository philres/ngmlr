#include <iostream>
#include <seqan/basic.h>
using namespace seqan;

int main()
{
	Dna a = 'a';
	std::cout << a << std::endl; 

	Dna5 b = 'f'; 
	std::cout << b << std::endl; 

	b = a;
	std::cout << b << std::endl; 

	Iupac c = b;
	std::cout << c << std::endl; 

	return 0;
}
