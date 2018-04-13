/**
    @file
    @author Rodrigo R. Paim
    @date 17/08/2013
    @section DESCRIPTION
    Header file of a simple graph implementation. For more information, look at Graph description.
*/

#ifndef GRAPH_H
#define GRAPH_H

#include "edge.h"

/**
This class is responsible for storing a graph with #_nodes vertices having IDs between 0 and n − 1 and its _nEdges labels. It can be directed (even though it is not fully implemented in some parts) or not, labeled or not. Moreover, each one of the _nLabels labels is a string, although it is not necessarily the best strategy due to high memory consumption. Each string is only stored once in #_labels to amortize this problem, even though multiple nodes may have it as a label. If the graph is undirected, each edge is stored twice.


Edges are stored sequentially in a unique giant array #__edges__ to optimize memory requests and try to have more cache hits. On the other hand, all the edge manipulations are done with #_edges.
*/
class Graph
{
public:
    int * __edges__; /**< Contiguous memory storage of the edges. Should not be referenced directly */

	std::map<std::string , int> _nodeIDtoint;
    std::map<int , std::string> _inttonodeID;
    int             _nodes; /**< Number of nodes in the graph */
    int *           _nEdges; /**< Degree of each vertex */
    int **          _edges; /**< Array of pointers to the edges of each node. Each pointer points to the first position the edges of a vertex is stored in #__edges__. Use always this instead of #_edges. For example, if you want to pick the first edge of node 2, you have to call #_edges[ 1 ][ 0 ] */

protected:
    bool            _directed; /**<Boolean variable to store whether the graph is directed (true) or not.*/


    /**
        Reserves memory for a graph with @a n nodes. It basically allocates memory for #_nEdges and #_edges
        @param n Value for the number of vertices
    */
	void    initialize( int n );
	std::string createMatchID( int i, int j );
	std::string createPairID( int i, int j, int k, int l);
	void untokenizeMatchID( std::string str, int & first, int & second );

public:
    /**
        Default constructor. Sets all int parameters to 0 and points all the arrays to NULL
    */
    //std::map< int, std::set<int> > getOverlappingCommunities();



    Graph( const Graph & g );

     Graph( bool directed = true);


    /**
        Destructor: Cleans #_labels and #_nodeLabel if #_nLabels is away from 0; cleans #__edges__, #_edges and #_nEdges if _nodes is away from 0. If these previous variables are 0, then it does not do anything to avoid Segmentation Fault
    */
    virtual ~Graph(); //Destructor


    /**
        Getter for #_nodes
        @return #_nodes
    */
    int             getNNodes();
	std::map<std::string , int> getnodeIDtoint();
    std::map<int , std::string> getinttonodeID();

    /**
        Calculates the number of edges in the graph
        @return @f$ \sum_{i=0}^{\_nodes-1} \_nEdges[ i ] @f$ if #_directed, @f$ \frac{1}{2}\sum_{i=0}^{\_nodes-1} \_nEdges[ i ] @f$ otherwise
    */
    int             getNEdges();

    /**
        Getter for #_labeled
        @return #_labeled
    */

    /**
        Getter for #_directed
        @return #_directed
    */
    bool            isDirected();

    /**
        Getter for the number of neighbors of a vertex
        @param index Target vertex
        @return #_nEdges[ index ]
    */
    int             getNNeighbors( int index );

    /**
        Getter for the array of neighbors of a vertex
        @param index Target vertex
        @return #_edges[ index ]
    */
    int *           getNeighbors( int index );

    /**
        Getter for the label of a vertex
        @param index Target vertex
        @return #_nodeLabel[ index ]
    */

    /**
        Method to check if an edge connecting @a u and @a v exists or not. As the edge set of a vertex is sorted, a binary search is performed, reducing the complexity to @f$log \; \delta @f$, where @f$\delta @f$ is the degree of @a u
        @param u Source vertex
        @param v Destination vertex
        @return #_edges[ u ][ binary_search( u, v ) ] == v
    */
    bool            edgeExists( int u, int v );

    /**
        Creates the graph from a list of Edge. This list must not contain Edge repetitions, because no check is performed. Iterates over all the entries once to discover the degree of each node and them store everything in the second iteration.
        @param nodes Number of nodes
        @param le List of Edge with the edges to be included in the graph
    */
    void            createGraph( int nodes, std::list< Edge > & le );

    /**
        Reads the graph edges from a text file. This file must be in a specific format: the first line has one int corresponding to #_nodes and all the other lines have two values, corresponding to the endpoints of the edges. These values must be in the range @f$[0,\_nodes-1]@f$ and edges cannot be repeated. If the graph is not directed, each edge will be stored twice, so the order is not important in this case (take care about this).
        @param filename Text file in the forementioned format with the graph structure
        @deprecated I/O cost makes this approach unfeasible for big graphs
    */
    void            readGraph( std::string filename );

    void            writeGraph( std::string output_file );




    /**
        Compares the edge set @f$ E @f$ of the object with the edge set @f$ E' @f$ of another graph and returns the Sørensen–Dice coefficient of them two. The graph to be compared with must have the same vertex set of the current object and each vertex must have the same unique ID in them both.
        @param g Graph whose edge set will be compared with the current object
        @return @f$ 2\frac{\vert E \cap E' \vert}{\vert E \vert + \vert E' \vert} @f$
    */
    double          compareGraph( const Graph & g );
    //void communityDetection( int display_level , bool do_renumber, bool verbose, double precision, std::string ocmfile);
    std::string retrieveString( char* buf );
    void display_time(const char *str);
    //void writeCommunities(std::string filename);
    //std::map<int, std::multiset<COMM, CompareCOMM> > seedCommunities(int nbNodes, std::list< Match> & seeds, std::map< int, int> leftComm, std::map< int, int> rightComm);
    //std::map<int, std::multiset<COMM, CompareCOMM> > updateSeedCommunities(std::map<int, std::multiset<COMM, CompareCOMM> > communities, int nbNodes, Match m, std::map< int, int> leftComm, std::map< int, int> rightComm);
    std::set<Match, CompareMatches> shmatikov(const Graph & g, std::list< Match> & seed, std::map< int, int> leftComm, std::map< int, int> rightComm,double ps);
    int  measure(int nodeindex, int level);
    void write_adj(std::string fname);

    std::set<Match, CompareMatches> align(const Graph & g, std::list< Match> & seed,int r, bool flag, int push_level, int push_size,std::string fnamei);
	std::set<Match, CompareMatches> proper(const Graph & g, std::list< Match> & seed,int r, bool tie_break);
};




#endif


