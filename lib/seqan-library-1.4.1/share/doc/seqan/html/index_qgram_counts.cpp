#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main ()
{
	StringSet< String<char> > mySet;
	resize(mySet, 4);
	mySet[0] = "tobeornottobe";
	mySet[1] = "thebeeonthecomb";
	mySet[2] = "hellobebe";
	mySet[3] = "beingjohnmalkovich";

	typedef Index< StringSet<String<char> >, IndexQGram<UngappedShape<2> > > TIndex;
	typedef Infix<Fibre<TIndex, QGramCounts>::Type const>::Type TCounts;

    TIndex myIndex(mySet);

    std::cout << "Number of sequences: " << countSequences(myIndex) << std::endl;  
    hash(indexShape(myIndex), "be");
    TCounts cnts = countOccurrencesMultiple(myIndex, indexShape(myIndex));
    for (unsigned i = 0; i < length(cnts); ++i)
       std::cout << cnts[i].i2 << " occurrences in sequence " << cnts[i].i1  << std::endl;

    
    String<double> distMat;
    getKmerSimilarityMatrix(myIndex,distMat);
    
    for( unsigned i=0; i < length(distMat); ++i)
        std::cout << distMat[i] << " ";
    std::cout << std::endl;
    
    return 0;
}
