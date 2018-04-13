#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>

#include "rg.h"

using namespace std;


void Erdos_Renyi::write_edges( int k, int n_structures, std::list< Edge > & le ) //O(M * logM)
{
    set< int > clique;
    set< pair< int, int > > edges;

    for ( int i = 0; i < n_structures; i++ )
    {
        clique.clear();
        while ( clique.size() < k )
        {
            clique.insert( rand() % _nodes );
        }

        for ( set< int >::iterator it = clique.begin(); it != clique.end(); ++it )
        {
            set< int >::iterator it2 = it; ++it2;
            for ( ; it2 != clique.end(); ++it2 )
            {
                edges.insert( make_pair( *it, *it2 ) );
            }
        }
    }

    for ( set< pair< int, int > >::iterator it = edges.begin(); it != edges.end(); ++it )
    {
        Edge e; e.u = it->first; e.v = it->second;
        le.push_back( e );
    }

    edges.clear();
}


Erdos_Renyi::Erdos_Renyi( int n ) : Graph( false )
{
    this->_nodes = n;
}

Erdos_Renyi::Erdos_Renyi( const Erdos_Renyi & er ) : Graph( ( const Graph & ) er )
{
}

void Erdos_Renyi::build( int avg, int k ) //O(M * logM)
{
    std::list< Edge > le;

    int n_structures = ( _nodes * avg ) / ( ( k-1 ) * k );
    write_edges( k, n_structures, le );
    createGraph( _nodes, le );
    le.clear();
}





Watts_Strogatz::Watts_Strogatz( int n ) : Graph( false)
{
    this->_nodes = n;
}

Watts_Strogatz::Watts_Strogatz( const Watts_Strogatz & ws ) : Graph( ( const Graph & ) ws )
{
}

void Watts_Strogatz::write_edges( int k, double beta, std::list< Edge > & le ) //k must be an even number
{
    set< Edges > edges;

    for ( int u = 0; u < _nodes; u++ ) //Wiring the lattice
    {
        for ( int j = 1; j <= k; j++ )
        {
            double p = double( rand() ) / RAND_MAX;
            int v;

            if ( beta > p ) //Change the edge
            {
                v = rand() % _nodes;
            }
            else //Add to a neighbor
            {
                v = ( u+j ) % _nodes;
            }

            while ( true )
            {
                while ( u == v ) //No self-loops
                {
                    v = rand() % _nodes;
                }

                int previous_size = edges.size();
                if ( u < v ) edges.insert( Edges(u, v ) );
                else edges.insert( Edges(v, u ) );
                if ( edges.size() > previous_size ) break; //New edge included!

                v = rand() % _nodes; //Choose another one
            }
        }
    }

    for ( set< Edges >::iterator it = edges.begin(); it != edges.end(); ++it ) //Adding to the return list
    {
        Edge e; e.u = it->getFrom(); e.v = it->getTo();
        le.push_back( e );
    }
}


void Watts_Strogatz::build( int k, double beta )
{
    std::list< Edge > le;
    write_edges( k, beta, le );
    createGraph( _nodes, le );
    le.clear();
}




Barabasi_Albert::Barabasi_Albert( int n ) : Graph( false)
{
    this->_nodes = n;
}

Barabasi_Albert::Barabasi_Albert( const Barabasi_Albert & ba ) : Graph( ( const Graph & ) ba )
{
}

void Barabasi_Albert::write_edges( int init, std::list< Edge > & le ) //O( M * logN + N² )
{
    int total_edges = 2 * init; //Number of initial edges (we count the two endpoints)
    double * linking_prob = new double[ _nodes ]; //Probability of attachment
    _nEdges = new int[ _nodes ]; //Number of edges of each node (used to compute the probability of attachment)

    //Create the network for the init+1 initial nodes
    for ( int i = 0; i < init; i++ )
    {
        _nEdges[ i ] = 1;
        Edge e; e.u = i; e.v = init;
        le.push_back( e );
    }
    _nEdges[ init ] = init;

    for ( int i = init+1; i < _nodes; i++ )
    {
        set< int > neighbors; //Neighbors of current node
        int sum = 0; //Cumulative sum of edges
        _nEdges[ i ] = 0;
        for ( int j = 0; j < i; j++ )
        {
            linking_prob[ j ] = double( sum ) / total_edges; //Changing the probability of attachment of the node
            sum += _nEdges[ j ];
        }

        int w = 0;
        while ( w < init ) //Connect the arriving node to init other ones
        {
            double p = double( rand() ) / RAND_MAX; //Random number
            int v = binary_search( i, linking_prob, p ); //Discover the node v that should be connected to i
            if ( neighbors.find( v ) != neighbors.end() ) continue; //Edge already exists

            neighbors.insert( v ); //Add v to the neighborhood so as not to avoid the creation of the same edge in the future
            _nEdges[ i ]++;
            _nEdges[ v ]++;
            total_edges += 2;
            Edge e; e.u = v; e.v = i;
            le.push_back( e );
            w++;
        }
    }

    delete [] linking_prob;
    delete [] _nEdges;

    _nEdges = NULL;
}

void Barabasi_Albert::build( int init ) //O( M * logN + N² )
{
    std::list< Edge > le;
    write_edges( init, le );
    createGraph( _nodes, le );
    le.clear();
}
