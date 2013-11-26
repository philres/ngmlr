//#include <iostream>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/align.h>

int main()
{
	using namespace seqan;
	typedef Value<Gaps<Dna5String, ArrayGaps> >::Type TValue;
	using namespace seqan;

    typedef String<Dna> TSequence;
    TSequence seq1 = "atcgaatgcgga";
    TSequence seq2 = "actcgttgca";
    Score<int> scoringScheme(0, -1, -1, -2);
    Align<TSequence, ArrayGaps> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);

    int score = globalAlignment(align, scoringScheme);
    std::cout << "Score = " << score << std::endl;
    std::cout << align << std::endl;
    score = globalAlignment(align, MyersHirschberg());
    std::cout << "Score = " << score << std::endl;
    std::cout << align << std::endl;
    typedef StringSet<TSequence, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

    TStringSet string_set;
    appendValue(string_set, seq1);
    appendValue(string_set, seq2);
    TAlignmentGraph alignment_graph(string_set);

    score = globalAlignment(alignment_graph, scoringScheme, Gotoh());
    std::cout << "Score = " << score << std::endl;
    std::cout << alignment_graph << std::endl;
    return 0;
}
