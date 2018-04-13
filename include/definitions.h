#pragma once

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <list>
#include <map>
#include <string>
#include <utility>

#define SEPARATOR "$" /**< String separator used as a token for tokenization and untokenization */

#define TMP_DIR "./_tmp" /**< @deprecated Used for file operations. Should be removed in next version! */
#define TMP_EGO TMP_DIR "/tmp_ego" /**< @deprecated Used for file operations. Should be removed in next version! */
#define TMP_LABEL TMP_DIR "/tmp_label" /**< @deprecated Used for file operations. Should be removed in next version! */
#define TMP_GRAPH TMP_DIR "/tmp_graph" /**< @deprecated Used for file operations. Should be removed in next version! */

/** Exceptions that can be raised in the middle of the code. Normally, whenever one of them occurs, the filename and line are displayed */
enum EXCEPTIONS {
                    NO_GRAPH, /**< Graph not created when method was called */
                    ALLOC_EXCEPTION, /**< Memory reallocation */
                    SEGFAULT /**< Common Segmentation Fault */
                };
/**
    Simple structure to store an edge where the endpoints have integer values.
*/
struct Edge
{
    int u, /**< Source endpoint */
        v; /**< Destination endpoint */
};


struct Match
{
	int lnode;
	int rnode;
	int value;
};

/*struct Tuple
{
	std::vector<int> vec;
	int value;
};

struct CompareTuples
{
  bool operator()(const Tuple & lhs, const Tuple & rhs) {
	if(lhs.value != rhs.value)
	{
		return lhs.value > rhs.value;
	}
	else
	{
		for(int i = 0; i < lhs.vec.size(); i++)
		{
			
		}
	}
		
  }
};*/



struct BlastVal
{
    std::string s; 
	int sp1;
	int sp2;
    int val;
};

struct CompareBlastVals
{
	bool operator()(const BlastVal & lhs, const BlastVal & rhs) {
		 if(lhs.val == rhs.val)
		{
			return lhs.s > rhs.s;
		}
		else
		{
			return lhs.val > rhs.val;
		}
	}

};






struct CompareMatches
{
  bool operator()(const Match & lhs, const Match & rhs) {
	  if(lhs.value == rhs.value)
	  {
		  if(lhs.lnode == rhs.lnode)
		  {
			  return ((((lhs.rnode * 0xf7f7f7f7) ^ 0x8364abf7) * 0xf00bf00b) ^ 0xf81bc437) > ((((rhs.rnode * 0xf7f7f7f7) ^ 0x8364abf7) * 0xf00bf00b) ^ 0xf81bc437);
		  }
		  return ((((lhs.lnode ^ 0xf7f7f7f7) * 0x8364abf7) ^ 0xf00bf00b) * 0xf81bc437) < ((((rhs.lnode ^ 0xf7f7f7f7) * 0x8364abf7) ^ 0xf00bf00b) * 0xf81bc437 );
	  }
	  else
	  {
		  return lhs.value > rhs.value;
	  }
  }
};

struct ProtPair
{
    int prot1,
		sp1,
        prot2,
		sp2,
		value; 
};



struct CompareProtPairs
{
  bool operator()(const ProtPair & lhs, const ProtPair & rhs) {
	  if(lhs.value == rhs.value)
	  {
		  if(lhs.prot1 == rhs.prot1)
		  {
			if(lhs.prot2 == rhs.prot2)
			{
				if(lhs.sp1 == rhs.sp1)
				{
					return (lhs.sp2 < rhs.sp2);
				}
				return (lhs.sp1 < rhs.sp1);
			}
			return (lhs.prot2 > rhs.prot2);
		  }
		  return (lhs.prot1 > rhs.prot1);
	  }
	  else
	  {
		  return (lhs.value > rhs.value);
	  }
  }
};


struct COMM
{
	int comm;
	int value;
};


struct CompareCOMM
{
bool operator()(const COMM & lhs, const COMM & rhs) {
	  if(lhs.value == rhs.value)
	  {
		  return (((((lhs.comm) ^ 0xf7f7f7f7) * 0x8364abf7) ^ 0xf00bf00b) * 0xf81bc437) < (((((rhs.comm) ^ 0xf7f7f7f7) * 0x8364abf7) ^ 0xf00bf00b) * 0xf81bc437 );;
	  }
	  else
	  {
		  return lhs.value > rhs.value;
	  }
  }
};
struct NodeVal
{
	int rnode;
	int value;
};
struct CompareNodeVal
{
bool operator()(const NodeVal & lhs, const NodeVal & rhs) {
	  if(lhs.value == rhs.value)
	  {
		  return lhs.rnode > rhs.rnode;
	  }
	  else
	  {
		  return lhs.value > rhs.value;
	  }
  }
};
/**
    Template function to perform a binary search in a sorted array.
    @param arr_size Length of @a arr
    @param arr Array where the search is performed. It must be sorted in descending order!
    @param val Value to be searched
    @return First index @a idx of @a arr such that @f$ arr[ idx ] \leq val @f$
*/
template <typename T>
int binary_search( int arr_size, T * arr, T val )
{
    int start, middle, end;

    start   = 0;
    end     = arr_size-1;
    middle  = ( start + end )/2;

    if ( !( arr[ start ] < val ) )
    {
        return start;
    }

    if ( !( val < arr[ end ] ) )
    {
        return end;
    }

    while ( start < end )
    {
        if ( arr[ middle ] < val )
        {
            start = middle;
        }
        else if ( val < arr[ middle ] )
        {
            end = middle;
        }
        else
        {
            return middle;
        }

        if ( middle == ( start + end )/2 )
        {
            break;
        }

        middle = ( start + end )/2;
    }

    return middle;
}

/**
    Template function to perform the union of two multiset in the form of list. Each one of the input lists must be sorted in descending order on the key value. If an element has cardinality @a 2 in one list, it must appear @a twice in the same list, for example.
    @param arr_1 First multiset
    @param arr_2 Second multiset
    @return Number of elements in the union of the two @a arr_1 and @a arr_2
*/
template <typename T>
int multiset_union( const std::list< T > & arr_1, const std::list< T > & arr_2 )
{
    int ret = 0;
    typename std::list< T >::const_iterator it_1 = arr_1.begin(), it_2 = arr_2.begin();

    while ( it_1 != arr_1.end() or it_2 != arr_2.end() )
    {
        ret++;

        if ( it_1 == arr_1.end() ) ++it_2;
        else if ( it_2 == arr_2.end() ) ++it_1;
        else if ( * it_1 < * it_2 ) ++it_1;
        else if ( * it_2 < * it_1 ) ++it_2;
        else { ++it_1; ++it_2; }
    }
    return ret;
}


/**
    Template function to perform the intersection of two multiset in the form of list. Each one of the input lists must be sorted in descending order on the key value. If an element has cardinality @a 2 in one list, it must appear @a twice in the same list, for example.
    @param arr_1 First multiset
    @param arr_2 Second multiset
    @return Number of elements in the intersection of the two @a arr_1 and @a arr_2
*/
template <typename T>
int multiset_intersection( const std::list< T > & arr_1, const std::list< T > & arr_2 )
{
    int ret = 0;
    typename std::list< T >::const_iterator it_1 = arr_1.begin(), it_2 = arr_2.begin();

    while ( it_1 != arr_1.end() or it_2 != arr_2.end() )
    {
        if ( it_1 == arr_1.end() or it_2 == arr_2.end() ) break;
        else if ( * it_1 < * it_2 ) ++it_1;
        else if ( * it_2 < * it_1 ) ++it_2;
        else { ++it_1; ++it_2; ret++; }
    }
    return ret;
}

#endif

