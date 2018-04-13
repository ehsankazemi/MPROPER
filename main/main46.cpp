#include "rg.h"
#include "alignment.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <ctime>

enum Type { ErdosRenyi, BarabasiAlbert, WattsStrogatz }; //0 for GNP, 1 for BA and 2 for WS

using namespace std;

bool matchCompareUniqueList(Match lhs, Match rhs) {
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


bool matchCompareSortList(Match lhs, Match rhs) {
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


int main( int argc, char * argv[] )
{
	srand( time( NULL ) );
	string sp1 = argv[1];
	string graph_file1 = "proc/" + sp1 + ".txt";
	string sp2 = argv[2];
	string graph_file2 = "proc/" + sp2 + ".txt";
	string sp3 = argv[3];
	string graph_file3 = "proc/" + sp3 + ".txt";
	string sp4 = argv[4];
	string graph_file4 = "proc/" + sp4 + ".txt";
	string sp5 = argv[5];
	string graph_file5 = "proc/" + sp5 + ".txt";

	

	int min_blast = atoi( argv[6] );
	int nb_seed = atoi( argv[7] );
	int int_flag = atoi( argv[8] );
	bool flag;
	if(int_flag == 1)
		flag = true;
	else
		flag = false;
	string filename = "scripts/match-" + std::string(argv[8]) + std::string("-") + sp1 + "-" + sp2 + "-" + sp3 + "-" + sp4 + "-" + sp5 + "-" + argv[6] + ".txt";
	ofstream output( filename.c_str() );
	cout << filename << endl;
	string blastfile = "scripts/hh_"+ sp1 + "_" +sp2 + "_" + sp3 + "_" + sp4 + "_" + sp5  + "_" + std::string(argv[6]) + ".txt";
	int order;
	//return 0;
	///////////////////////////////////////////////////////////////////////////////////////////
	//read graphs
	Graph * g1 = new Graph( false);
	g1->readGraph( graph_file1 );
	cout << g1->getNNodes() << "\t" << g1->getNEdges() << endl;
	
	Graph * g2 = new Graph( false);
	g2->readGraph( graph_file2 );
	cout << g2->getNNodes() << "\t" << g2->getNEdges() << endl;

	Graph * g3 = new Graph( false);
	g3->readGraph( graph_file3 );
	cout << g3->getNNodes() << "\t" << g3->getNEdges() << endl;

	Graph * g4 = new Graph( false);
	g4->readGraph( graph_file4 );
	cout << g4->getNNodes() << "\t" << g4->getNEdges() << endl;

	Graph * g5 = new Graph( false);
	g5->readGraph( graph_file5 );
	cout << g5->getNNodes() << "\t" << g5->getNEdges() << endl;


	Alignment * align = new Alignment();
	///////////////////////////////////////////////////////////////////////////////////////////
	//read seeds
	std::list< Quintuple > seed;
	std::set<Quintuple, CompareQuintuples> check;
	int u, v, w, x, y;
	float n;
	int node1, node2, node3, node4;
	ifstream file;
   	Quintuple Q;
	
	cout << blastfile << endl;
	file.open( blastfile.c_str() );
	int values[4];
	int cnt = 0;
	while (file >> u >> v >> w >> x >> y)
	{	
		int n = 200;
		if(cnt > nb_seed)
			break;
		cnt++;		
		Q.node1 = u; Q.node2 = v; Q.node3 = w; Q.node4 = x; Q.node5 = y; Q.value = (int)(n + 1);
		std::set<Quintuple, CompareQuintuples>::iterator it = check.find(Q);
		if(it != check.end())
		{
			cout << u << "\t" << y << "\t" << v << "\t" << x << "\t"  << w << "\t" << n << endl;
			cout << "++++++++++++++" << endl;
			return 0;
		}
		seed.push_back(Q);
		check.insert(Q);

	}
	file.close();
	cout << check.size() << " seed Done " << seed.size() << endl;
	///////////////////////////////////////////////////////////////////////////////////////////

	Graph * g = new Graph [5];
	g[0] = *g1;
	g[1] = *g2;
	g[2] = *g3;
	g[3] = *g4;
	g[4] = *g5;

	cout << "seed: " << seed.size() << endl;
	std::set<Quintuple, CompareQuintuples> matches5 = align->multiple(g, seed, flag);
	cout << matches5.size() << endl;
	for (std::set<Quintuple, CompareQuintuples>::iterator it = matches5.begin(); it != matches5.end(); ++it )
	{
		output << it->node1 << "\t" << it->node2 << "\t" << it->node3 << "\t" << it->node4 << "\t" << it->node5 << std::endl;
	}

	
	
	delete g1;
	delete g2;
	delete g3;
	delete g4;
	delete g5;
	seed.clear();
	matches5.clear();
	output.close();
	return 0;
}
