#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;


int main() {
    typedef Graph<Directed<> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    typedef Size<TGraph>::Type TSize;
    TSize numEdges = 10;
    TVertexDescriptor edges[] = {0,1, 0,4, 1,2, 1,4, 2,3, 2,4, 4,1, 4,5, 5,2, 5,3};
    TGraph g;
    addEdges(g,edges, numEdges);
    std::cout << g << std::endl;
    String<unsigned int> capMap;    
    unsigned int capacity[] =    {16,  13,  12,  10,  20,  9,   4,   14,  7,   4};
    assignEdgeMap(g,capMap, capacity);
    String<unsigned int> flow;
    unsigned int valF = fordFulkersonAlgorithm(g, 0, 3, capMap, flow);
    std::cout << "Ford-Fulkerson (Value of the flow = " << valF << ")" << std::endl;
    TEdgeIterator itEdge(g);
    for(;!atEnd(itEdge);goNext(itEdge)) {
        std::cout << "(" << sourceVertex(itEdge) << "," << targetVertex(itEdge) << "): ";
        std::cout << "Flow: " << getProperty(flow, getValue(itEdge)) << ", Capacity: "
                    << getProperty(capMap, getValue(itEdge)) << std::endl;
	}
	return 0;
}
