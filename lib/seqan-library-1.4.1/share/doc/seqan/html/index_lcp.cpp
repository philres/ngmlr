#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
	StringSet< String<char> > mySet;
	resize(mySet, 3);
	mySet[0] = "SeqAn is a library for sequence analysis.";
	mySet[1] = "The String class is the fundamental sequence type in SeqAn.";
	mySet[2] = "Subsequences can be handled with SeqAn's Segment class.";

	typedef Index< StringSet<String<char> > > TMyIndex;
	TMyIndex myIndex(mySet);

    indexRequire(myIndex, EsaLcp() );
    indexRequire(myIndex, EsaSA());
    
    for ( Size<TMyIndex>::Type i=0; i<length(myIndex); ++i){
        SAValue<TMyIndex>::Type p = saAt(i,myIndex);
        std::cout << i << " " << lcpAt(i,myIndex) << " " << p << " " << suffix(mySet,p) << std::endl;
    }

	return 0;
}
