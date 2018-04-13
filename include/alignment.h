#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "edge.h"
#include "graph.h"
using namespace std;

class Alignment
{
private:


public:
	Alignment();
	std::string createPairID( int i, int j, int k, int l);
	std::string createID( int i, int j);
	std::map <int , vector<int> > pair(Graph g[], std::list<Match> & seed);
	std::map <int, vector<int>> multipleK(Graph g[], int nb_networks, std::map <int , vector<int> > seed, bool check_flag);
	virtual ~Alignment(); //Destructor


};







#endif
