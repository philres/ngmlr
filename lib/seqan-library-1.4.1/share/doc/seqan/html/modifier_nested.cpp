#include <iostream>
#include <seqan/file.h>
#include <seqan/modifier.h>

using namespace seqan;

int main ()
{
	String<Dna> myString = "attacgg";
	typedef ModifiedString<String<Dna>, ModComplementDna>   TMyComplement;
	typedef ModifiedString<TMyComplement, ModReverse>       TMyReverseComplement;

	TMyReverseComplement myReverseComplement(myString);
	std::cout << myString << std::endl;
	std::cout << myReverseComplement << std::endl;
	replace(myString, 1, 1, "cgt");
	std::cout << myString << std::endl;
	std::cout << myReverseComplement << std::endl;
	std::cout << DnaStringReverseComplement(myString) << std::endl;
	return 0;
}
