#include <iostream>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    String<char> myString = "How many wood would a woodchuck chuck. A woodchuck chucks as much wood as a woodchuck could";

    typedef Index<String<char> > TMyIndex;
    TMyIndex myIndex(myString);

    Iterator<TMyIndex, TopDown<ParentLinks<PreorderEmptyEdges> > >::Type tdIterator(myIndex);
    Size<TMyIndex>::Type count;

    while (!atEnd(tdIterator))
    {
        count = countChildren(tdIterator);
        if (count >= 3)
        {
            std::cout << "Representative " << representative(tdIterator) << " has " <<  count << " children  and ";
            std::cout << countOccurrences(tdIterator) << " occurrences " << std::endl;
        }
        if (isLeaf(tdIterator))
            std::cout << "The node is a leaf " << std::endl;

        tdIterator++;
    }

    return 0;
}
