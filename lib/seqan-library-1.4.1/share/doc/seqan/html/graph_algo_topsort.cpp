#include <iostream>
#include <seqan/graph_algorithms.h>


using namespace seqan;


int main() {
	typedef Graph<Directed<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;
	TSize numEdges = 9;
	TVertexDescriptor edges[] = {0,3, 0,1, 1,2, 3,2, 5,7, 5,6, 6,7, 6,3, 8,7};
	TGraph g;
	addEdges(g, edges, numEdges);
	std::cout << g << std::endl;
	String<std::string> nameMap;
	std::string names[] = {"shirt", "tie", "jacket", "belt", "watch", "undershorts", "pants", "shoes", "socks"};
	assignVertexMap(g,nameMap, names);
	String<TVertexDescriptor> order;
	topologicalSort(g, order);
	std::cout << "Topological sort: " << std::endl;
	typedef Iterator<String<TVertexDescriptor> >::Type TStringIterator;
	TStringIterator it = begin(order);
	TStringIterator itEnd = end(order);
	while(it != itEnd) {
		std::cout << getProperty(nameMap, getValue(it)) << ",";
		goNext(it);
	}
	std::cout << std::endl;
	return 0;
}
