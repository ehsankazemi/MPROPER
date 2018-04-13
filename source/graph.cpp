#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include "graph.h"
using namespace std;
#define EdgeSeparator SEPARATOR /**< Token used to separate strings and ints */


bool matchCompareUnique(Match lhs, Match rhs) {
	  if(lhs.value == rhs.value)
	  {
		  if(lhs.lnode == rhs.lnode)
		  {
			  return (lhs.rnode == rhs.rnode);
		  }
		  return (lhs.lnode == rhs.lnode);
	  }
	  else
	  {
		  return lhs.value == rhs.value;
	  }
  }


bool matchCompareSort(Match lhs, Match rhs) {
	  if(lhs.value == rhs.value)
	  {
		  if(lhs.lnode == rhs.lnode)
		  {
			  return (lhs.rnode > rhs.rnode);
		  }
		  return (lhs.lnode > rhs.lnode);
	  }
	  else
	  {
		  return lhs.value > rhs.value;
	  }
  }


Graph::Graph( bool directed)
{
    _nodes      = 0;
    _nEdges     = NULL;
    _edges      = NULL;
    __edges__   = NULL;

    _directed   = directed;
}

std::map<std::string , int> Graph::getnodeIDtoint()
{
	return _nodeIDtoint;
}

std::map<int , std::string> Graph::getinttonodeID()
{
	return _inttonodeID;
}

Graph::Graph( const Graph & g )
{

    _nodes     = g._nodes;
    _nEdges    = new int[ _nodes ];
    _edges     = new int * [ _nodes ];

    int total_edges = 0; //Total edges is used to create a contiguous block of memory

    for ( int i = 0; i < _nodes; i++ )
    {
        _nEdges[ i ]    = g._nEdges[ i ];
        total_edges    += _nEdges[ i ];
    }

    __edges__ = new int[ total_edges ];

    total_edges = 0;
    for ( int i = 0; i < _nodes; i++ )
    {
        this->_edges[ i ] = __edges__ + total_edges; //Setting the pointers to the list of edges to the right place
        for ( int j = 0; j < _nEdges[ i ]; j++ ) //Copying the edges
        {
            this->_edges[ i ][ j ] = g._edges[ i ][ j ];
        }
        total_edges += _nEdges[ i ];
    }

    _directed  = g._directed;

}





std::string Graph::createMatchID( int i, int j )
{
    std::stringstream sstm;
    sstm << i << "$" << j;
    return sstm.str();
}

std::string Graph::createPairID( int i, int j, int k, int l )
{
    std::stringstream sstm;
    sstm << i << "$" << j << "$" << k << "$" << l;
    return sstm.str();
}




void Graph::untokenizeMatchID( std::string str, int & first, int & second )
{
    std::string ES = EdgeSeparator;
    int index = str.find( ES );
    first  = atoi(str.substr( 0, index ).c_str());
    if ( index != -1 )
        second = atoi(str.substr( str.find( ES ) + ES.size() ).c_str());
}

void Graph::createGraph( int nodes, list< Edge > & le )
{
	initialize( nodes );
    for ( list< Edge >::iterator it = le.begin(); it != le.end(); ++it ) //Calculating the number of neighbors of each node
    {
        _nEdges[ it->u ]++;
        if ( !_directed )
        {
            _nEdges[ it->v ]++;
        }
    }

    //Create a contiguous memory buffer for storing the edges
    long long size = 0;
    for ( int i = 0; i < nodes; i++ )
    {
        size += _nEdges[ i ];
    }
    __edges__ = new int[ size ];
    size = 0;
    for ( int i = 0; i < nodes; i++ )
    {
        _edges[ i ] = __edges__ + size;
        size += _nEdges[ i ];
    }

    //Read for the second time: writing edges: No need for try-catch block. File exists unless user deleted it...
    int * pos = new int[ nodes ]; //How many edges have been added for each node
    for ( int i = 0; i < nodes; i++ )
    {
        pos[ i ] = 0;
    }

    for ( list< Edge >::iterator it = le.begin(); it != le.end(); ++it ) //Writing the edges
    {
        int u = it->u, v = it->v;
        _edges[ u ][ pos[ u ]++ ] = v;
        if ( ! _directed )
        {
            _edges[ v ][ pos[ v ]++ ] = u;
        }
    }

   for ( int i = 0; i < _nodes; i++ ) //Sort the neighbors to simplify search in the future
    {
        if ( _nEdges[ i ] > 0 )
        {
        	//random_shuffle( _edges[ i ], _edges[ i ] + _nEdges[ i ] );
        	sort( _edges[ i ], _edges[ i ] + _nEdges[ i ] );
        }
    }

    delete [] pos;
}

Graph::~Graph()
{
    if ( _nodes > 0 )
    {
        delete [] __edges__;
        delete [] _edges;
        delete [] _nEdges;
    }
}

void Graph::initialize( int n )
{
    _nodes     = n;
    _nEdges    = new int[ n ];
    _edges     = new int * [ n ];

    for ( int i = 0; i < n; i++ )
    {
        _nEdges[ i ]   = 0;
        _edges[ i ]    = NULL;
    }
    _directed   = false;
}


int Graph::getNEdges()
{
    int ret = 0;
    for ( int i = 0; i < _nodes; i++ )
    {
        ret += _nEdges[ i ];
    }
    if ( _directed ) //If it is not a directed graph, then it means that each edge was counted twice
    {
        return ret;
    }
    return ret/2;
}


void Graph::readGraph( string filename )
{
    int u, v, n;
    string e1, e2;

    if ( __edges__ != 0 ) //If the graph has been already constructed, avoid a memory leak
    {
        cout << __FILE__ << " " << __LINE__ << endl;
        throw ALLOC_EXCEPTION;
    }

    ifstream file;
    ifstream file_orig;
    
    try
    {
        //Change the graph format 
        file_orig.open( filename.c_str() );
        int nb_nodes = 0;
        while ( file_orig >> e1 >> e2 )
        {
            if(_nodeIDtoint.find(e1) == _nodeIDtoint.end())
            {
                _nodeIDtoint[e1] = nb_nodes;
                _inttonodeID[nb_nodes] = e1;
                nb_nodes++;
            }
            if(_nodeIDtoint.find(e2) == _nodeIDtoint.end())
            {
                _nodeIDtoint[e2] = nb_nodes;
                _inttonodeID[nb_nodes] = e2;
                nb_nodes++;
            }
        }
        file_orig.close();
        //Read for the first time: discover number of edges
        file.open( filename.c_str() );
        n = nb_nodes;
        initialize( n ); //Initialize the graph with its correct size
        while ( file >> e1 >> e2 )
        {
            u = _nodeIDtoint[e1];
            v = _nodeIDtoint[e2];           
            _nEdges[ u ]++;
            if ( !_directed )
            {
                _nEdges[ v ]++;
            }
        }
        file.close();
    }
    catch(...)
    {
        cout << "Problem reading the graph file: " << __FILE__ << " " << __LINE__ << endl;
        return;
    }

    //Create a contiguous memory buffer for storing the edges
    int size = 0;
    for ( int i = 0; i < n; i++ )
    {
        size += _nEdges[ i ];
    }
    __edges__ = new int[ size ];
    size = 0;
    for ( int i = 0; i < n; i++ )
    {
        _edges[ i ] = __edges__ + size;
        size += _nEdges[ i ];
    }

    //Read for the second time: writing edges: No need for try-catch block. File exists unless user deleted it...
    int * it = new int[ n ]; //How many edges have been added for each node
    for ( int i = 0; i < n; i++ )
    {
        it[ i ] = 0;
    }

    file.open( filename.c_str() );
     while ( file >> e1 >> e2 )
      {
        u = _nodeIDtoint[e1];
        v = _nodeIDtoint[e2];  
        _edges[ u ][ it[ u ]++ ] = v;
        if ( !_directed ) //If it is not a directed graph, add the opposite edge too
        {
            _edges[ v ][ it[ v ]++ ] = u;
        }
    }
    file.close();

    for ( int i = 0; i < _nodes; i++ ) //Sort the neighbors to simplify search in the future
    {
        if ( _nEdges[ i ] > 0 )
            sort( _edges[ i ], _edges[ i ] + _nEdges[ i ] );
    }

    delete [] it;
}



int Graph::getNNodes()
{
    return _nodes;
}


int Graph::getNNeighbors( int index )
{
    if ( index >= _nodes )
    {
        return 0;
    }
    return _nEdges[ index ];
}


int* Graph::getNeighbors( int index )
{
    if ( index >= _nodes )
    {
        return NULL;
    }
    return _edges[ index ];
}




bool Graph::isDirected()
{
    return _directed;
}




bool Graph::edgeExists( int u, int v ) //O( logN )
{
    if ( _nEdges[ u ] == 0 or _nEdges[ v ] == 0 ) return false;
    int pos = binary_search( _nEdges[ u ], _edges[ u ], v );
    if ( pos < 0 or pos >= _nodes )
    {
        cout << __FILE__ << " " << __LINE__ << endl;
        throw SEGFAULT;
    }
    return _edges[ u ][ pos ] == v;
}


void Graph::writeGraph( string output_file )
{
    ofstream output( output_file.c_str() );
    output << _nodes << endl; //Output number of nodes
    for ( int i = 0; i < _nodes; i++ )
    {
        for ( int j = 0; j < _nEdges[ i ]; j++ )
        {
            if ( i > _edges[ i ][ j ] ) //Prints each edge only once
                continue;
            output << i << " " << _edges[ i ][ j ] << endl;
        }
    }
    output.close();
}

double Graph::compareGraph( const Graph & g ) //O( M )
{
    //Compare two graphs according to the existence (or not) of an edge in both graphs
    int total_size = 0, total_match = 0;
    int * it_this, * it_g;

    for ( int u = 0; u < _nodes; u++ ) //O( M )
    {
        it_this = _edges[ u ];
        it_g    = g._edges[ u ];

        total_size += ( _nEdges[ u ] + g._nEdges[ u ] );

        //The neighbors of both graphs are sorted, so we can perform a matching in linear time
        while ( it_this != _edges[ u ] + _nEdges[ u ] and it_g != g._edges[ u ] + g._nEdges[ u ] )
        {
            if ( *it_this < *it_g )
            {
                it_this++;
            }
            else if ( *it_this > *it_g )
            {
                it_g++;
            }
            else
            {
                it_this++; it_g++;
                total_match += 2;
            }
        }
    }

    if ( total_size == 0 ) return 1.; //If there is no edge in both graphs, then it means that they match exactly, as they have the same number of neighbors

    return double( total_match ) / double( total_size ) ;
}

void Graph::display_time(const char *str) {
  time_t rawtime;
  time ( &rawtime );
  cerr << str << ": " << ctime (&rawtime);
}

std::string Graph::retrieveString( char* buf ) {

    size_t len = 0;
    while (buf[ len ] != '\0')
    {
        len++;
    }
    return string( buf, len );
}


