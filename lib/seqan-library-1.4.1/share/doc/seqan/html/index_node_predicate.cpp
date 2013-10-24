#include <iostream>
#include <seqan/index.h>

using namespace seqan;

struct TMyConstraints {
	double p_min;
	unsigned int replen_max;
	bool _cachedPred;
};

namespace seqan 
{
	// custom TSpec for our customized wotd-Index
	struct TMyIndex;

	template <typename TText>
	struct Cargo<Index<TText, IndexWotd<TMyIndex> > > {
		typedef TMyConstraints Type;
	};

	// node predicate
	template <typename TText, typename TSpec>
	bool nodePredicate(Iter<Index<TText, IndexWotd<TMyIndex> >, TSpec> &it) 
	{
		TMyConstraints &cons = cargo(container(it));
		unsigned int delta = countSequences(container(it)) * repLength(it);
		unsigned int textLen = length(container(it));

		if (textLen <= delta) return cons._cachedPred = true;
		return cons._cachedPred = 
			((double)countOccurrences(it) / (double)(textLen - delta)) > cons.p_min;
	}

	// monotonic hull
	template <typename TText, typename TSpec>
	bool nodeHullPredicate(Iter<Index<TText, IndexWotd<TMyIndex> >, TSpec> &it) 
	{
		TMyConstraints const &cons = cargo(container(it));
		unsigned int textLen = length(container(it));

		if (repLength(it) > cons.replen_max) return false;

		unsigned int delta = countSequences(container(it)) * cons.replen_max;
		if (textLen <= delta) return true;
		return ((double)countOccurrences(it) / 
				(double)(textLen - delta)) > cons.p_min;
	}
}

int main ()
{
	String<char> myString = "How many wood would a woodchuck chuck.";

	typedef Index< String<char>, IndexWotd<TMyIndex> > TMyIndex;
	TMyIndex myIndex(myString);

	cargo(myIndex).replen_max = 10;
	cargo(myIndex).p_min = 0.05;

	typedef Iterator< TMyIndex, TopDown<ParentLinks<> > >::Type TConstrainedIterator;
	TConstrainedIterator myConstrainedIterator(myIndex);

	goBegin(myConstrainedIterator);
	while (!atEnd(myConstrainedIterator))
	{

		std::cout << countOccurrences(myConstrainedIterator) << "x  ";

		std::cout << "\t\"" << representative(myConstrainedIterator) << '\"' << std::endl;

		goNext(myConstrainedIterator);
	}

	return 0;
}
