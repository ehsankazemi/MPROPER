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
#include "alignment.h"
using namespace std;
Alignment::Alignment(void)
{
    cout << "Object is being created" << endl;
}

Alignment::~Alignment()
{
    matches.clear();
}


std::string Alignment::createPairID( int i, int j, int k, int l )
{
    std::stringstream sstm;
	if(j < l)
	    sstm << i << "$" << j << "$" << k << "$" << l;
	else
		sstm << k << "$" << l << "$" << i << "$" << j;
    return sstm.str();
}

std::string Alignment::createID( int i, int j)
{
    std::stringstream sstm;
    sstm << i << "$" << j;
    return sstm.str();
}


std::map <int, vector<int>> Alignment::pair(Graph g[], std::list<Match> & seed)
{
	std::set<Match, CompareMatches> matches;
	std::map< std::string, int> score_map;
	std::multiset<ProtPair, CompareProtPairs> scores;
	int nb_networks = 2;
	std::set<int> * matched_nodes = new std::set<int> [nb_networks];//set of matched nodes
	int check = 0;
	int t_cnt = 0;
	std::map< int, std::set<int>>* matched_neis =new std::map< int, std::set<int>> [nb_networks]; //matched neighbors of a node
	std::map< int, int>** alignments =new std::map< int, int>*[nb_networks]; //pairwise alignments
		for (int i = 0; i < nb_networks; ++i)
	{
	   alignments[i] = new std::map< int, int> [nb_networks]; //pairwise alignments
	}
	int cnt = 0;
	int tuple[nb_networks];
	int nb_tuple = 0;
	std::map <int , vector<int> > tuples;
	std::map<int, int> * node_tuple = new std::map<int, int> [nb_networks];

	for ( list< Match>::iterator it = seed.begin(); it != seed.end(); ++it )
	{	
		Match q; q.lnode = it->lnode; q.rnode = it->rnode; q.value = 2;
		for(int i = 0; i < nb_networks; i++)
			tuple[i] = -1;
		tuple[0] = it->lnode; tuple[1] = it->rnode;
		for(int i = 0; i < nb_networks; i++){
			for(int j = i+1; j < nb_networks; j++){
				if(tuple[i] != -1 and tuple[j] != -1){
					alignments[i][j][tuple[i]] = tuple[j];
					alignments[j][i][tuple[j]] = tuple[i];
				}
			}
		}
		for(int sp_nb = 0; sp_nb < nb_networks; sp_nb++)
		{			
			if(tuple[sp_nb] != -1)
			{
				for(int i=0; i < g[sp_nb]._nEdges[ tuple[sp_nb] ]; i++)
				{		
					matched_neis[sp_nb][g[sp_nb]._edges[tuple[sp_nb]][i]].insert(tuple[sp_nb]); //add tuple[sp_nb] to all the matched_neis of its neighbors
				}
			}
		}
		matched_nodes[0].insert(it->lnode);
		matched_nodes[1].insert(it->rnode);
		std::vector<int> v;
		for(int i = 0; i < nb_networks; i++)
		{
			v.push_back(tuple[i]);
			if(tuple[i] != -1)
			{
				node_tuple[i][tuple[i]] = nb_tuple;
			}	
		}
		tuples[nb_tuple] = v;
		nb_tuple++;
		cnt += 1;
	}
	//return tuples;
	//////////////////////////////////////////////////////////////////////////
	for (list<Match>::iterator it = seed.begin(); it != seed.end(); ++it )
	{
		int node1 = it->lnode;int node2 = it->rnode;
		std::list<int> nodeneis[nb_networks];
		int node_array[nb_networks]; node_array[0] = node1; node_array[1] = node2;
		for(int sp1 = 0; sp1 < nb_networks; sp1++)
		{
			for(int sp2 = sp1+1; sp2 < nb_networks; sp2++)
			{
				if(node_array[sp1] != -1 and node_array[sp2] != -1)
				{
					//cout << node_array[sp1] << "\t" << node_array[sp2] << endl;
					for(int i=0; i < g[sp1]._nEdges[node_array[sp1]]; i++)
					{
						for(int j=0; j < g[sp2]._nEdges[node_array[sp2]]; j++)
						{
						if(matched_nodes[sp1].find(g[sp1]._edges[node_array[sp1]][i]) == matched_nodes[sp1].end() && matched_nodes[sp2].find(g[sp2]._edges[node_array[sp2]][j]) == matched_nodes[sp2].end())
						{
							int node1 = g[sp1]._edges[node_array[sp1]][i];
							int node2 = g[sp2]._edges[node_array[sp2]][j];
							std::string ID = createPairID(node1, sp1, node2, sp2);
							if(score_map.find(ID) == score_map.end())
							{
								score_map[ID] = 1;
								ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 1;
								scores.insert(pp);
							}
							else
							{
								int val = score_map[ID];
								ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = val;
								std::multiset<ProtPair, CompareProtPairs>::iterator it = scores.find(pp);
								if(it != scores.end())
									scores.erase(it);
								val+=1;
								score_map[ID] = val;
								pp.value = val;
								scores.insert(pp);
							}	
						}
						}
					}
				}
			}

		}
	}
	
	/////////////////////////////////////////////////////////////////////////
	
	std::multiset<ProtPair, CompareProtPairs>::iterator it = scores.begin();
	//cout << scores.size() << endl;
	while(scores.size() > 0 and it->value > 1)
	{	
		//cout << it->sp1 << "\t" << it->sp2 << "\t" << it->value << endl;
		int node1 = it->prot1;
		int node2 = it->prot2;
		int sp1 = it->sp1;
		int sp2 = it->sp2;
		bool spread_flag = false;
		std::vector<ProtPair> spread_pairs;
		if(node_tuple[sp1].find(node1) == node_tuple[sp1].end() and node_tuple[sp2].find(node2) == node_tuple[sp2].end())
		{
			std::vector<int> v;
			for(int i = 0; i < nb_networks; i++)
			{
				v.push_back(-1);
			}
			if(node1 == -1 or node2 == -1)
			{
				cout << "hey" << endl;
			}
			v[sp1] = node1;
			v[sp2] = node2;
			node_tuple[sp1][node1] = nb_tuple;
			node_tuple[sp2][node2] = nb_tuple;
			tuples[nb_tuple] = v;
			matched_nodes[sp1].insert(node1);
			matched_nodes[sp2].insert(node2);
			nb_tuple++;
			spread_flag = true;
			ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 0;
			spread_pairs.push_back(pp);
		}

		if(node_tuple[sp1].find(node1) != node_tuple[sp1].end() and node_tuple[sp2].find(node2) == node_tuple[sp2].end())
		{
			int tuple_ind = node_tuple[sp1][node1];
			if(tuples[tuple_ind][sp2] == -1)
			{
				std::vector<int> v;
				v = tuples[tuple_ind];
				v[sp2] = node2;
				tuples[tuple_ind] = v;
				node_tuple[sp2][node2] = tuple_ind;
				matched_nodes[sp2].insert(node2);
				spread_flag = true;
				for(int i = 0; i < nb_networks; i++)
				{
					if(i != sp2 and tuples[tuple_ind][i] != -1)
					{
						ProtPair pp; pp.prot1 = tuples[tuple_ind][i]; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 0;
						spread_pairs.push_back(pp);
					}
				}
			}
		}
		
		if(node_tuple[sp1].find(node1) == node_tuple[sp1].end() and node_tuple[sp2].find(node2) != node_tuple[sp2].end())
		{
			int tuple_ind = node_tuple[sp2][node2];
			if(tuples[tuple_ind][sp1] == -1)
			{
				std::vector<int> v;
				v = tuples[tuple_ind];
				v[sp1] = node1;
				tuples[tuple_ind] = v;
				node_tuple[sp1][node1] = tuple_ind;
				matched_nodes[sp1].insert(node1);
				//spread_flag = true;
				
			}
		}
		if(node_tuple[sp1].find(node1) != node_tuple[sp1].end() and node_tuple[sp2].find(node2) != node_tuple[sp2].end())
		{
			int tuple_ind1 = node_tuple[sp1][node1];
			int tuple_ind2 = node_tuple[sp2][node2];
			if(tuples[tuple_ind1][sp2] == -1 and tuples[tuple_ind2][sp1] == -1)
			{		
				if(tuple_ind1 < tuple_ind2 or true)
				{
					for(int i = 0; i < nb_networks; i++)
					{
						if(tuples[tuple_ind1][i] == -1 and tuples[tuple_ind2][i] != -1)
						{
							tuples[tuple_ind1][i] = tuples[tuple_ind2][i];
							tuples[tuple_ind2][i] = -1;
							node_tuple[i][tuples[tuple_ind1][i]] = tuple_ind1;
						}
					}
				}
			}
			else
			{
				int cnt1 = 0;
				int cnt2 = 0;
				int cntt = 0;
				for(int i = 0; i < nb_networks; i++)
				{
					if(tuples[tuple_ind1][i] != -1) cnt1++;
					if(tuples[tuple_ind2][i] != -1) cnt2++;
					if(tuples[tuple_ind1][i] != -1 or tuples[tuple_ind2][i] != -1) cntt++;
				}
			}
		}
		if(spread_flag)
		{		
			for(std::vector<ProtPair>::iterator itv = spread_pairs.begin(); itv != spread_pairs.end(); ++itv)
    		{
				int node1 = itv->prot1;
				int node2 = itv->prot2;
				int sp1 = itv->sp1;
				int sp2 = itv->sp2;
				for(int i=0; i < g[sp1]._nEdges[node1]; i++)
				{
					for(int j=0; j < g[sp2]._nEdges[node2]; j++)
					{
						if(matched_nodes[sp1].find(g[sp1]._edges[node1][i]) == matched_nodes[sp1].end() && matched_nodes[sp2].find(g[sp2]._edges[node2][j]) == matched_nodes[sp2].end())
						{
					
							int tnode1 = g[sp1]._edges[node1][i];
							int tnode2 = g[sp2]._edges[node2][j];
							std::string ID = createPairID(tnode1, sp1, tnode2, sp2);
							if(score_map.find(ID) == score_map.end())
							{						
								score_map[ID] = 1;
								ProtPair pp; pp.prot1 = tnode1; pp.sp1 = sp1; pp.prot2 = tnode2; pp.sp2 = sp2; pp.value = 1;
								scores.insert(pp);
							}
							else
							{							
								int val = score_map[ID];
								ProtPair pp; pp.prot1 = tnode1; pp.sp1 = sp1; pp.prot2 = tnode2; pp.sp2 = sp2; pp.value = val;
								std::multiset<ProtPair, CompareProtPairs>::iterator itt = scores.find(pp);
								if(itt != scores.end())
									scores.erase(itt);
								val+=1;
								score_map[ID] = val;
								pp.value = val;
								scores.insert(pp);
							}	
						}
					}
				}
			}
			it = scores.begin();
		}
		else
		{
			it = scores.erase(it);
			std::string ID = createPairID(node1, sp1, node2, sp2);
			score_map.erase(ID);
		}
	}
	//cout << "size: " << tuples.size() << endl;
	return tuples;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

std::set<Quintuple, CompareQuintuples> Alignment::multiple(Graph g[], std::list<Quintuple> & seed, bool flag)
{

	std::set<Quintuple, CompareQuintuples> matches;
	std::list< Quintuple > seedtmp;
	std::map< std::string, int> score_map;
	std::multiset<ProtPair, CompareProtPairs> scores;
	std::set<int> * matched_nodes = new std::set<int> [5];//set of matched nodes
	std::set<std::string> seedcheck;
	int check = 0;
	int t_cnt = 0;
	std::map< int, std::set<int>>* matched_neis =new std::map< int, std::set<int>> [5]; //matched neighbors of a node
	std::map< int, int>** alignments =new std::map< int, int>*[5]; //pairwise alignments
	for (int i = 0; i < 5; ++i)
	{
	   alignments[i] = new std::map< int, int> [5]; //pairwise alignments
	}
	int cnt = 0;
	int tuple[5];
	int nb_tuple = 0;
	std::map <int , vector<int> > tuples;
	std::map<int, int> * node_tuple = new std::map<int, int> [5];

	for ( list< Quintuple >::iterator it = seed.begin(); it != seed.end(); ++it )
	{	
		Quintuple q; q.node1 = it->node1; q.node2 = it->node2; q.node3 = it->node3; q.node4 = it->node4; q.node5 = it->node5; q.value = 6;
		for(int i = 0; i < 5; i++)
			tuple[i] = -1;
		tuple[0] = it->node1;tuple[1] = it->node2;tuple[2] = it->node3;tuple[3] = it->node4;tuple[4] = it->node5;
		for(int i = 0; i < 5; i++){
			for(int j = i+1; j < 5; j++){
				if(tuple[i] != -1 and tuple[j] != -1){
					alignments[i][j][tuple[i]] = tuple[j];
					alignments[j][i][tuple[j]] = tuple[i];
				}
			}
		}
		//cout << "here!" << endl;
		for(int sp_nb = 0; sp_nb < 5; sp_nb++)
		{
			//cout << "here! " << sp_nb << endl;			
			if(tuple[sp_nb] != -1)
			{
				for(int i=0; i < g[sp_nb]._nEdges[ tuple[sp_nb] ]; i++)
				{		
					matched_neis[sp_nb][g[sp_nb]._edges[tuple[sp_nb]][i]].insert(tuple[sp_nb]); //add tuple[sp_nb] to all the matched_neis of its neighbors
				}
			}
		}
		//cout << "out!" << endl;
		matched_nodes[0].insert(it->node1);
		matched_nodes[1].insert(it->node2);
		matched_nodes[2].insert(it->node3);
		matched_nodes[3].insert(it->node4);
		matched_nodes[4].insert(it->node5);
		std::vector<int> v;
		for(int i = 0; i < 5; i++)
		{
			v.push_back(tuple[i]);
			if(tuple[i] != -1)
			{
				node_tuple[i][tuple[i]] = nb_tuple;
			}	
		}
		tuples[nb_tuple] = v;
		//int tuple_ind = nb_tuple;
		//cout << tuples[tuple_ind][0] << "\t" << tuples[tuple_ind][1] << "\t" << tuples[tuple_ind][2] << "\t" << tuples[tuple_ind][3] << "\t" << tuples[tuple_ind][4] << endl;
		nb_tuple++;
		//matches.insert(q);//add all the seeds to the final tuples
		cnt += 1;
	}
	if(not flag)
		return matches;

	


	for (list<Quintuple>::iterator it = seed.begin(); it != seed.end(); ++it )
	{
		int node1 = it->node1;int node2 = it->node2;int node3 = it->node3;int node4 = it->node4;int node5 = it->node5;
		std::list<int> nodeneis[5];
		int node_array[5]; node_array[0] = node1; node_array[1] = node2; node_array[2] = node3; node_array[3] = node4; node_array[4] = node5;
		for(int sp1 = 0; sp1 < 5; sp1++)
		{
			for(int sp2 = sp1+1; sp2 < 5; sp2++)
			{
				if(node_array[sp1] != -1 and node_array[sp2] != -1)
				{
					//cout << node_array[sp1] << "\t" << node_array[sp2] << endl;
					for(int i=0; i < g[sp1]._nEdges[node_array[sp1]]; i++)
					{
						for(int j=0; j < g[sp2]._nEdges[node_array[sp2]]; j++)
						{
						if(matched_nodes[sp1].find(g[sp1]._edges[node_array[sp1]][i]) == matched_nodes[sp1].end() && matched_nodes[sp2].find(g[sp2]._edges[node_array[sp2]][j]) == matched_nodes[sp2].end())
						{
							int node1 = g[sp1]._edges[node_array[sp1]][i];
							int node2 = g[sp2]._edges[node_array[sp2]][j];
							std::string ID = createPairID(node1, sp1, node2, sp2);
							if(score_map.find(ID) == score_map.end())
							{
								score_map[ID] = 1;
								ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 1;
								scores.insert(pp);
							}
							else
							{
								int val = score_map[ID];
								ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = val;
								std::multiset<ProtPair, CompareProtPairs>::iterator it = scores.find(pp);
								if(it != scores.end())
									scores.erase(it);
								val+=1;
								score_map[ID] = val;
								pp.value = val;
								scores.insert(pp);
							}	
						}
						}
					}
				}
			}

		}

	}
	std::multiset<ProtPair, CompareProtPairs>::iterator it = scores.begin();
	//std::map <int , vector<int> > tuples;
	//std::map<int, int> * node_tuple = new std::map<int, int> [5];
	//while(it != scores.end() and it->value > 0)
	while(scores.size() > 0)
	{	
		//cout << it->sp1 << "\t" << it->sp2 << "\t" << it->value << endl;
		int node1 = it->prot1;
		int node2 = it->prot2;
		int sp1 = it->sp1;
		int sp2 = it->sp2;
		bool spread_flag = false;
		if(node_tuple[sp1].find(node1) == node_tuple[sp1].end() and node_tuple[sp2].find(node2) == node_tuple[sp2].end())
		{
			std::vector<int> v;
			for(int i = 0; i < 5; i++)
			{
				v.push_back(-1);
			}
			if(node1 == -1 or node2 == -1)
			{
				cout << "hey" << endl;
			}
			v[sp1] = node1;
			v[sp2] = node2;
			//cout << v[0] << "\t" << v[1] << "\t" << v[2] << "\t" << v[3] << "\t" << v[4] << endl;
			node_tuple[sp1][node1] = nb_tuple;
			node_tuple[sp2][node2] = nb_tuple;
			tuples[nb_tuple] = v;
			matched_nodes[sp1].insert(node1);
			matched_nodes[sp2].insert(node2);
			nb_tuple++;
			spread_flag = true;
			//cout << node1 << "\t" << sp1 << endl;
		}

		if(node_tuple[sp1].find(node1) != node_tuple[sp1].end() and node_tuple[sp2].find(node2) == node_tuple[sp2].end())
		{
			int tuple_ind = node_tuple[sp1][node1];
			if(tuples[tuple_ind][sp2] == -1)
			{
				//cout << tuples[tuple_ind][0] << "\t" << tuples[tuple_ind][1] << "\t" << tuples[tuple_ind][2] << "\t" << tuples[tuple_ind][3] << "\t" << tuples[tuple_ind][4] << endl;
				std::vector<int> v;
				v = tuples[tuple_ind];
				v[sp2] = node2;
				tuples[tuple_ind] = v;
				//cout << tuples[tuple_ind][0] << "\t" << tuples[tuple_ind][1] << "\t" << tuples[tuple_ind][2] << "\t" << tuples[tuple_ind][3] << "\t" << tuples[tuple_ind][4] << endl;	
				//cout << "+++++" << endl;
				node_tuple[sp2][node2] = tuple_ind;
				matched_nodes[sp2].insert(node2);
				spread_flag = true;
			}
		}
		
		if(node_tuple[sp1].find(node1) == node_tuple[sp1].end() and node_tuple[sp2].find(node2) != node_tuple[sp2].end())
		{
			int tuple_ind = node_tuple[sp2][node2];
			if(tuples[tuple_ind][sp1] == -1)
			{
				//cout << tuples[tuple_ind][0] << "\t" << tuples[tuple_ind][1] << "\t" << tuples[tuple_ind][2] << "\t" << tuples[tuple_ind][3] << "\t" << tuples[tuple_ind][4] << endl;
				std::vector<int> v;
				v = tuples[tuple_ind];
				v[sp1] = node1;
				tuples[tuple_ind] = v;
				//cout << tuples[tuple_ind][0] << "\t" << tuples[tuple_ind][1] << "\t" << tuples[tuple_ind][2] << "\t" << tuples[tuple_ind][3] << "\t" << tuples[tuple_ind][4] << endl;	
				//cout << "-----" << endl;
				node_tuple[sp1][node1] = tuple_ind;
				matched_nodes[sp1].insert(node1);
				spread_flag = true;
			}
		}
		if(node_tuple[sp1].find(node1) != node_tuple[sp1].end() and node_tuple[sp2].find(node2) != node_tuple[sp2].end())
		{
			int tuple_ind1 = node_tuple[sp1][node1];
			int tuple_ind2 = node_tuple[sp2][node2];
			if(tuples[tuple_ind1][sp2] == -1 and tuples[tuple_ind2][sp1] == -1)
			{
				bool merge_flag = true;
				for(int i = 0; i < 5; i++)
				{
					if(tuples[tuple_ind1][i] != -1 and tuples[tuple_ind2][i] != -1)
						merge_flag = false;
				}
				if(merge_flag or true)
				{		
					for(int i = 0; i < 5; i++)
					{
						if(tuples[tuple_ind1][i] == -1 and tuples[tuple_ind2][i] != -1)
						{
							tuples[tuple_ind1][i] = tuples[tuple_ind2][i];
							tuples[tuple_ind2][i] = -1;
							node_tuple[i][tuples[tuple_ind1][i]] = tuple_ind1;
						}
					}
				}
			}
		}
		if(spread_flag)
		{
			//cout << "spread!" << endl;			
			for(int i=0; i < g[sp1]._nEdges[node1]; i++)
			{
				for(int j=0; j < g[sp2]._nEdges[node2]; j++)
				{
					if(matched_nodes[sp1].find(g[sp1]._edges[node1][i]) == matched_nodes[sp1].end() && matched_nodes[sp2].find(g[sp2]._edges[node2][j]) == matched_nodes[sp2].end())
					{
						//cout << "check!" << endl;						
						int tnode1 = g[sp1]._edges[node1][i];
						int tnode2 = g[sp2]._edges[node2][j];
						std::string ID = createPairID(tnode1, sp1, tnode2, sp2);
						if(score_map.find(ID) == score_map.end())
						{
							//cout << "else!" << endl;							
							score_map[ID] = 1;
							ProtPair pp; pp.prot1 = tnode1; pp.sp1 = sp1; pp.prot2 = tnode2; pp.sp2 = sp2; pp.value = 1;
							scores.insert(pp);
						}
						else
						{
							//cout << "else!" << endl;								
							int val = score_map[ID];
							ProtPair pp; pp.prot1 = tnode1; pp.sp1 = sp1; pp.prot2 = tnode2; pp.sp2 = sp2; pp.value = val;
							std::multiset<ProtPair, CompareProtPairs>::iterator itt = scores.find(pp);
							if(itt != scores.end())
								scores.erase(itt);
							val+=1;
							score_map[ID] = val;
							pp.value = val;
							scores.insert(pp);
						}	
					}
				}
			}
			it = scores.begin();
			//cout << "start!" << endl;
		}
		else
		{
			it = scores.erase(it);
			std::string ID = createPairID(node1, sp1, node2, sp2);
			score_map.erase(ID);
		}
	}
	cout << "size " << tuples.size() << endl;
	//return matches;
	int five_cnt = 0;
	for(int i = 0; i < tuples.size(); i++)
	{
		int cnt = 0;
		for(int j = 0; j < 5; j++)
		{
			if(tuples[i][j] != -1)
			{
				cnt++;
			}
		}
		/*if(cnt == 5)
		{
			five_cnt++;		
			cout << five_cnt << "\t" << tuples[i][0] << "\t" << tuples[i][1] << "\t" << tuples[i][2] << "\t" << tuples[i][3] << "\t" << tuples[i][4] << endl;
		}*/
		Quintuple q; q.node1 = tuples[i][0]; q.node2 = tuples[i][1]; q.node3 = tuples[i][2]; q.node4 = tuples[i][3]; q.node5 = tuples[i][4]; q.value = cnt;
		if(cnt > 1)
		{
			matches.insert(q);
		}
	}
	return matches;
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
std::map <int, vector<int>> Alignment::five(Graph g[], std::list<Quintuple> & seed, bool seed_flag, bool check_flag)
{
	cout << "start" << endl;
	std::map< std::string, int> score_map;
	std::multiset<ProtPair, CompareProtPairs> scores;
	std::set<std::string> seedcheck;
	std::set<std::string> matched;
	int nb_networks = 5;
	std::set<int> * matched_nodes = new std::set<int> [nb_networks];//set of matched nodes
	int check = 0;
	int t_cnt = 0;
	int cnt = 0;
	int tuple[nb_networks];
	int nb_tuple = 0;
	std::map <int , vector<int> > tuples;
	
	std::map<int, int> * node_tuple = new std::map<int, int> [nb_networks];
	cout << "spread0 "  << seed.size() << endl;
	for ( list<Quintuple>::iterator it = seed.begin(); it != seed.end(); ++it )
	{	
		Quintuple q; q.node1 = it->node1; q.node2 = it->node2; q.node3 = it->node3; q.node4 = it->node4; q.node5 = it->node5; q.value = 5;
		for(int i = 0; i < nb_networks; i++)
			tuple[i] = -1;
		tuple[0] = it->node1; tuple[1] = it->node2; tuple[2] = it->node3; tuple[3] = it->node4; tuple[4] = it->node5;

		matched_nodes[0].insert(it->node1);
		std::string ID = createID(it->node1, 0);
		matched.insert(ID);
		matched_nodes[1].insert(it->node2);
		ID = createID(it->node2, 1);
		matched.insert(ID);
		matched_nodes[2].insert(it->node3);
		ID = createID(it->node3, 2);
		matched.insert(ID);
		matched_nodes[3].insert(it->node4);
		ID = createID(it->node4, 3);
		matched.insert(ID);
		matched_nodes[4].insert(it->node5);
		ID = createID(it->node5, 4);
		matched.insert(ID);
		std::vector<int> v;
		for(int i = 0; i < nb_networks; i++)
		{
			v.push_back(tuple[i]);
			if(tuple[i] != -1)
			{
				node_tuple[i][tuple[i]] = nb_tuple;
			}	
		}
		tuples[nb_tuple] = v;
		nb_tuple++;
		cnt += 1;
	}
	if(not seed_flag)
		return tuples;
	cout << "spread1" << endl;
	//////////////////////////////////////////////////////////////////////////
	for (list<Quintuple>::iterator it = seed.begin(); it != seed.end(); ++it )
	{
		int node1 = it->node1; int node2 = it->node2; int node3 = it->node3; int node4 = it->node4; int node5 = it->node5;
		int node_array[nb_networks]; node_array[0] = node1; node_array[1] = node2;  node_array[2] = node3; node_array[3] = node4; node_array[4] = node5;
		for(int sp1 = 0; sp1 < nb_networks; sp1++)
		{
			for(int sp2 = sp1+1; sp2 < nb_networks; sp2++)
			{				
				if(node_array[sp1] != -1 and node_array[sp2] != -1)
				{
					std::string IDD = createPairID(node_array[sp1], sp1, node_array[sp2], sp2);
					if(seedcheck.find(IDD) == seedcheck.end())
					{
						seedcheck.insert(IDD);						
						for(int i=0; i < g[sp1]._nEdges[node_array[sp1]]; i++)
						{
							for(int j=0; j < g[sp2]._nEdges[node_array[sp2]]; j++)
							{
								std::string ID1 = createID(g[sp1]._edges[node_array[sp1]][i], sp1);
								std::string ID2 = createID(g[sp2]._edges[node_array[sp2]][j], sp2);
								if(check_flag || (matched.find(ID1) == matched.end() && matched.find(ID2) == matched.end()))
								{
									int node1 = g[sp1]._edges[node_array[sp1]][i];
									int node2 = g[sp2]._edges[node_array[sp2]][j];
									std::string ID = createPairID(node1, sp1, node2, sp2);
									if(score_map.find(ID) == score_map.end())
									{
										score_map[ID] = 1;
										ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 1;
										scores.insert(pp);
									}
									else
									{
										int val = score_map[ID];
										ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = val;
										std::multiset<ProtPair, CompareProtPairs>::iterator it = scores.find(pp);
										if(it != scores.end())
											scores.erase(it);
										val+=1;
										score_map[ID] = val;
										pp.value = val;
										scores.insert(pp);
									}	
								}
							}
						}
					}
				}
			}

		}
	}
	cout << "define2" << endl;
	/////////////////////////////////////////////////////////////////////////
	
	std::multiset<ProtPair, CompareProtPairs>::iterator it = scores.begin();
	cout << "score size: " << scores.size() << endl;
	while(scores.size() > 0 and it->value > 0)
	{	
		
		int node1 = it->prot1;
		int node2 = it->prot2;
		int sp1 = it->sp1;
		int sp2 = it->sp2;
		
		bool spread_flag = false;
		std::vector<ProtPair> spread_pairs;
		if(node_tuple[sp1].find(node1) == node_tuple[sp1].end() and node_tuple[sp2].find(node2) == node_tuple[sp2].end())
		{
			//if(node1 == node2)
				//cout << seedcheck.size() << " ++++ " << node1 << "\t" << node2 << "\t" << sp1 << "\t" << sp2 << endl;;			
			std::vector<int> v;
			for(int i = 0; i < nb_networks; i++)
			{
				v.push_back(-1);
			}
			if(node1 == -1 or node2 == -1)
			{
				cout << "hey" << endl;
			}
			v[sp1] = node1;
			v[sp2] = node2;
			node_tuple[sp1][node1] = nb_tuple;
			node_tuple[sp2][node2] = nb_tuple;
			tuples[nb_tuple] = v;
			//matched_nodes[sp1].insert(node1);
			//matched_nodes[sp2].insert(node2);
			//std::string ID = createID(node1, sp1);
			//matched.insert(ID);
			//createID(node2, sp2);
			//matched.insert(ID);
			nb_tuple++;
			spread_flag = true;
			ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 1;
			spread_pairs.push_back(pp);
		}

		if(node_tuple[sp1].find(node1) != node_tuple[sp1].end() and node_tuple[sp2].find(node2) == node_tuple[sp2].end())
		{
			int tuple_ind = node_tuple[sp1][node1];
			if(tuples[tuple_ind][sp2] == -1)
			{
				std::vector<int> v;
				v = tuples[tuple_ind];
				v[sp2] = node2;
				tuples[tuple_ind] = v;
				node_tuple[sp2][node2] = tuple_ind;
				matched_nodes[sp2].insert(node2);
				//std::string ID = createID(node2, sp2);
				//matched.insert(ID);
				spread_flag = true;
				////////////////////////////////////////////////////////////////////
				for(int i = 0; i < nb_networks; i++)
				{
					if(i != sp2 and tuples[tuple_ind][i] != -1)
					{
						ProtPair pp; pp.prot1 = tuples[tuple_ind][i]; pp.sp1 = i; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 0;
						spread_pairs.push_back(pp);
					}
				}
				/*ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 0;
				spread_pairs.push_back(pp);*/
				////////////////////////////////////////////////////////////////////
			}
		}
		
		if(node_tuple[sp1].find(node1) == node_tuple[sp1].end() and node_tuple[sp2].find(node2) != node_tuple[sp2].end())
		{
			int tuple_ind = node_tuple[sp2][node2];
			if(tuples[tuple_ind][sp1] == -1)
			{
				std::vector<int> v;
				v = tuples[tuple_ind];
				v[sp1] = node1;
				tuples[tuple_ind] = v;
				node_tuple[sp1][node1] = tuple_ind;
				matched_nodes[sp1].insert(node1);
				//std::string ID = createID(node1, sp1);
				//matched.insert(ID);
				spread_flag = true;
				///////////////////////////////////////////////////////////////////////
				for(int i = 0; i < nb_networks; i++)
				{
					if(i != sp1 and tuples[tuple_ind][i] != -1)
					{
						ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = tuples[tuple_ind][i]; pp.sp2 = i; pp.value = 0;
						spread_pairs.push_back(pp);
					}
				}
				/*ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 0;
				spread_pairs.push_back(pp);*/
				//////////////////////////////////////////////////////////////////////
			}
		}
	
		if(node_tuple[sp1].find(node1) != node_tuple[sp1].end() and node_tuple[sp2].find(node2) != node_tuple[sp2].end())
		{
			int tuple_ind1 = node_tuple[sp1][node1];
			int tuple_ind2 = node_tuple[sp2][node2];
			if(tuples[tuple_ind1][sp2] == -1 and tuples[tuple_ind2][sp1] == -1)
			{		
				bool merge_flag = true;
				for(int i = 0; i < 5; i++)
				{
					if(tuples[tuple_ind1][i] != -1 and tuples[tuple_ind2][i] != -1)
						merge_flag = false;
				}
				if(merge_flag)
				{
					spread_flag = true;					
					for(int i = 0; i < nb_networks; i++)
					{
						for(int j = 0; j < nb_networks; j++)
						{
							if(tuples[tuple_ind1][i] != -1 and tuples[tuple_ind2][j] != -1)
							{
								ProtPair pp; pp.prot1 = tuples[tuple_ind1][i]; pp.sp1 =i; pp.prot2 = tuples[tuple_ind2][j]; pp.sp2 = j; pp.value = 0;
								spread_pairs.push_back(pp);
							}
						}
					}
					for(int i = 0; i < 5; i++)
					{
						if(tuples[tuple_ind1][i] == -1 and tuples[tuple_ind2][i] != -1)
						{
							tuples[tuple_ind1][i] = tuples[tuple_ind2][i];
							tuples[tuple_ind2][i] = -1;
							node_tuple[i][tuples[tuple_ind1][i]] = tuple_ind1;
						}
					}
				}
			}
		}
		if(spread_flag)
		{		
			for(std::vector<ProtPair>::iterator itv = spread_pairs.begin(); itv != spread_pairs.end(); ++itv)
    		{			
				int node1 = itv->prot1;
				int node2 = itv->prot2;
				int sp1 = itv->sp1;
				int sp2 = itv->sp2;
				std::string IDD = createPairID(node1, sp1, node2, sp2);
				if(seedcheck.find(IDD) == seedcheck.end())
				{
					seedcheck.insert(IDD);
					for(int i=0; i < g[sp1]._nEdges[node1]; i++)
					{
						for(int j=0; j < g[sp2]._nEdges[node2]; j++)
						{
							int tnode1 = g[sp1]._edges[node1][i];
							int tnode2 = g[sp2]._edges[node2][j];							
							std::string ID1 = createID(tnode1, sp1);
							std::string ID2 = createID(tnode2, sp2);
							if(check_flag ||  (matched.find(ID1) == matched.end() && matched.find(ID2) == matched.end()))
							{		
								std::string ID = createPairID(tnode1, sp1, tnode2, sp2);
								if(score_map.find(ID) == score_map.end())
								{						
									score_map[ID] = 1;
									ProtPair pp; pp.prot1 = tnode1; pp.sp1 = sp1; pp.prot2 = tnode2; pp.sp2 = sp2; pp.value = 1;
									scores.insert(pp);
								}
								else
								{							
									int val = score_map[ID];
									ProtPair pp; pp.prot1 = tnode1; pp.sp1 = sp1; pp.prot2 = tnode2; pp.sp2 = sp2; pp.value = val;
									std::multiset<ProtPair, CompareProtPairs>::iterator itt = scores.find(pp);
									if(itt != scores.end())
										scores.erase(itt);
									val+=1;
									score_map[ID] = val;
									pp.value = val;
									scores.insert(pp);
								}	
							}
						}
					}
				}
			}
			it = scores.begin();
		}
		else
		{
			it = scores.erase(it);
			std::string ID = createPairID(node1, sp1, node2, sp2);
			score_map.erase(ID);
		}
	}
	cout << "size: " << tuples.size() << endl;
	return tuples;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::map <int, vector<int>> Alignment::multipleK(Graph g[], int nb_networks, std::map <int , vector<int> > seed, bool check_flag)
{
	cout << "start" << endl;
	std::map< std::string, int> score_map;
	std::multiset<ProtPair, CompareProtPairs> scores;
	std::set<std::string> seedcheck;
	std::set<std::string> matched;
	std::set<int> * matched_nodes = new std::set<int> [nb_networks];//set of matched nodes
	int check = 0;
	int t_cnt = 0;
	int cnt = 0;
	int nb_tuple = 0;
	std::map <int , vector<int> > tuples;
	
	std::map<int, int> * node_tuple = new std::map<int, int> [nb_networks];
	cout << "spread0 "  << seed.size() << endl;
	for(int seed_ind = 0; seed_ind <  seed.size(); seed_ind++)
	{	
		std::vector<int> tuple;		
		tuple = seed[seed_ind];
		for(int i = 0; i < nb_networks; i++)
		{
			if(tuple[i] != -1)
			{
				std::string ID = createID(tuple[i], 0);
				matched.insert(ID);
			}
		}	
		for(int i = 0; i < nb_networks; i++)
		{
			if(tuple[i] != -1)
			{
				node_tuple[i][tuple[i]] = nb_tuple;
			}	
		}
		tuples[nb_tuple] = tuple;
		nb_tuple++;
		cnt += 1;
	}
	cout << "spread1" << endl;
	//////////////////////////////////////////////////////////////////////////
	for(int seed_ind = 0; seed_ind <  seed.size(); seed_ind++)
	{	
		std::vector<int> tuple;
		tuple = seed[seed_ind];	
		for(int sp1 = 0; sp1 < nb_networks; sp1++)
		{
			for(int sp2 = sp1+1; sp2 < nb_networks; sp2++)
			{				
				if(tuple[sp1] != -1 and tuple[sp2] != -1)
				{
					std::string IDD = createPairID(tuple[sp1], sp1, tuple[sp2], sp2);
					if(seedcheck.find(IDD) == seedcheck.end())
					{
						seedcheck.insert(IDD);						
						for(int i=0; i < g[sp1]._nEdges[tuple[sp1]]; i++)
						{
							for(int j=0; j < g[sp2]._nEdges[tuple[sp2]]; j++)
							{
								std::string ID1 = createID(g[sp1]._edges[tuple[sp1]][i], sp1);
								std::string ID2 = createID(g[sp2]._edges[tuple[sp2]][j], sp2);
								if(check_flag || (matched.find(ID1) == matched.end() && matched.find(ID2) == matched.end()))
								{
									int node1 = g[sp1]._edges[tuple[sp1]][i];
									int node2 = g[sp2]._edges[tuple[sp2]][j];
									std::string ID = createPairID(node1, sp1, node2, sp2);
									if(score_map.find(ID) == score_map.end())
									{
										score_map[ID] = 1;
										ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 1;
										scores.insert(pp);
									}
									else
									{
										int val = score_map[ID];
										ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = val;
										std::multiset<ProtPair, CompareProtPairs>::iterator it = scores.find(pp);
										if(it != scores.end())
											scores.erase(it);
										val+=1;
										score_map[ID] = val;
										pp.value = val;
										scores.insert(pp);
									}	
								}
							}
						}
					}
				}
			}

		}
	}
	cout << "define2" << endl;
	cout << "score size: " << scores.size() << endl;
	for(int run = 0; run < 1; run++)
	{
		std::multiset<ProtPair, CompareProtPairs>::iterator it = scores.begin();		
		cout << "run: " << run << endl;		
		while(scores.size() > 0 and it->value > 1)
		{	
		
			int node1 = it->prot1;
			int node2 = it->prot2;
			int sp1 = it->sp1;
			int sp2 = it->sp2;
		
			bool spread_flag = false;
			std::vector<ProtPair> spread_pairs;
			if(node_tuple[sp1].find(node1) == node_tuple[sp1].end() and node_tuple[sp2].find(node2) == node_tuple[sp2].end())
			{			
				std::vector<int> v;
				for(int i = 0; i < nb_networks; i++)
				{
					v.push_back(-1);
				}
				if(node1 == -1 or node2 == -1)
				{
					cout << "hey" << endl;
				}
				v[sp1] = node1;
				v[sp2] = node2;
				node_tuple[sp1][node1] = nb_tuple;
				node_tuple[sp2][node2] = nb_tuple;
				tuples[nb_tuple] = v;
				matched_nodes[sp1].insert(node1);
				matched_nodes[sp2].insert(node2);
				std::string ID = createID(node1, sp1);
				matched.insert(ID);
				createID(node2, sp2);
				matched.insert(ID);
				nb_tuple++;
				spread_flag = true;
				ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 1;
				spread_pairs.push_back(pp);
			}

			if(node_tuple[sp1].find(node1) != node_tuple[sp1].end() and node_tuple[sp2].find(node2) == node_tuple[sp2].end())
			{
				int tuple_ind = node_tuple[sp1][node1];
				if(tuples[tuple_ind][sp2] == -1)
				{
					std::vector<int> v;
					v = tuples[tuple_ind];
					v[sp2] = node2;
					tuples[tuple_ind] = v;
					node_tuple[sp2][node2] = tuple_ind;
					matched_nodes[sp2].insert(node2);
					std::string ID = createID(node2, sp2);
					matched.insert(ID);
					spread_flag = true;
					////////////////////////////////////////////////////////////////////
					for(int i = 0; i < nb_networks; i++)
					{
						if(i != sp2 and tuples[tuple_ind][i] != -1)
						{
							ProtPair pp; pp.prot1 = tuples[tuple_ind][i]; pp.sp1 = i; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 0;
							spread_pairs.push_back(pp);
						}
					}
					/*ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 0;
					spread_pairs.push_back(pp);*/
					////////////////////////////////////////////////////////////////////
				}
			}
		
			if(node_tuple[sp1].find(node1) == node_tuple[sp1].end() and node_tuple[sp2].find(node2) != node_tuple[sp2].end())
			{
				int tuple_ind = node_tuple[sp2][node2];
				if(tuples[tuple_ind][sp1] == -1)
				{
					std::vector<int> v;
					v = tuples[tuple_ind];
					v[sp1] = node1;
					tuples[tuple_ind] = v;
					node_tuple[sp1][node1] = tuple_ind;
					matched_nodes[sp1].insert(node1);
					std::string ID = createID(node1, sp1);
					matched.insert(ID);
					spread_flag = true;
					///////////////////////////////////////////////////////////////////////
					for(int i = 0; i < nb_networks; i++)
					{
						if(i != sp1 and tuples[tuple_ind][i] != -1)
						{
							ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = tuples[tuple_ind][i]; pp.sp2 = i; pp.value = 0;
							spread_pairs.push_back(pp);
						}
					}
					/*ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 0;
					spread_pairs.push_back(pp);*/
					//////////////////////////////////////////////////////////////////////
				}
			}
	
			if(node_tuple[sp1].find(node1) != node_tuple[sp1].end() and node_tuple[sp2].find(node2) != node_tuple[sp2].end())
			{
				int tuple_ind1 = node_tuple[sp1][node1];
				int tuple_ind2 = node_tuple[sp2][node2];
				if(tuples[tuple_ind1][sp2] == -1 and tuples[tuple_ind2][sp1] == -1)
				{		
					bool merge_flag = true;
					for(int i = 0; i < nb_networks; i++)
					{
						if(tuples[tuple_ind1][i] != -1 and tuples[tuple_ind2][i] != -1)
							merge_flag = false;
					}
					if(merge_flag)
					{
						spread_flag = true;					
						for(int i = 0; i < nb_networks; i++)
						{
							for(int j = 0; j < nb_networks; j++)
							{
								if(tuples[tuple_ind1][i] != -1 and tuples[tuple_ind2][j] != -1)
								{
									ProtPair pp; pp.prot1 = tuples[tuple_ind1][i]; pp.sp1 =i; pp.prot2 = tuples[tuple_ind2][j]; pp.sp2 = j; pp.value = 0;
									spread_pairs.push_back(pp);
								}
							}
						}
						for(int i = 0; i < nb_networks; i++)
						{
							if(tuples[tuple_ind1][i] == -1 and tuples[tuple_ind2][i] != -1)
							{
								tuples[tuple_ind1][i] = tuples[tuple_ind2][i];
								tuples[tuple_ind2][i] = -1;
								node_tuple[i][tuples[tuple_ind1][i]] = tuple_ind1;
							}
						}
					}
				}
			}
			if(spread_flag)
			{		
				for(std::vector<ProtPair>::iterator itv = spread_pairs.begin(); itv != spread_pairs.end(); ++itv)
				{			
					int node1 = itv->prot1;
					int node2 = itv->prot2;
					int sp1 = itv->sp1;
					int sp2 = itv->sp2;
					std::string IDD = createPairID(node1, sp1, node2, sp2);
					if(seedcheck.find(IDD) == seedcheck.end())
					{
						seedcheck.insert(IDD);
						for(int i=0; i < g[sp1]._nEdges[node1]; i++)
						{
							for(int j=0; j < g[sp2]._nEdges[node2]; j++)
							{
								int tnode1 = g[sp1]._edges[node1][i];
								int tnode2 = g[sp2]._edges[node2][j];							
								std::string ID1 = createID(tnode1, sp1);
								std::string ID2 = createID(tnode2, sp2);
								if(check_flag ||  (matched.find(ID1) == matched.end() && matched.find(ID2) == matched.end()))
								{		
									std::string ID = createPairID(tnode1, sp1, tnode2, sp2);
									if(score_map.find(ID) == score_map.end())
									{						
										score_map[ID] = 1;
										ProtPair pp; pp.prot1 = tnode1; pp.sp1 = sp1; pp.prot2 = tnode2; pp.sp2 = sp2; pp.value = 1;
										scores.insert(pp);
									}
									else
									{							
										int val = score_map[ID];
										ProtPair pp; pp.prot1 = tnode1; pp.sp1 = sp1; pp.prot2 = tnode2; pp.sp2 = sp2; pp.value = val;
										std::multiset<ProtPair, CompareProtPairs>::iterator itt = scores.find(pp);
										if(itt != scores.end())
											scores.erase(itt);
										val+=1;
										score_map[ID] = val;
										pp.value = val;
										scores.insert(pp);
									}	
								}
							}
						}
					}
				}
				it = scores.begin();
			}
			else
			{
				it = scores.erase(it);
				std::string ID = createPairID(node1, sp1, node2, sp2);
				score_map.erase(ID);
			}
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/*if(run == 0)
		{
			std::multiset<ProtPair, CompareProtPairs>::iterator ite = scores.begin();
			std::vector<ProtPair> spread_pairs;
			while(ite != scores.end())
			{	
			
				int node1 = ite->prot1;
				int node2 = ite->prot2;
				int sp1 = ite->sp1;
				int sp2 = ite->sp2;
				ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 1;
				spread_pairs.push_back(pp);
				ite++;
			}
			int cnt_it = 0;
			for(std::vector<ProtPair>::iterator itv = spread_pairs.begin(); itv != spread_pairs.end(); ++itv)
			{			
				int node1 = itv->prot1;
				int node2 = itv->prot2;
				int sp1 = itv->sp1;
				int sp2 = itv->sp2;
				//cout << cnt_it << " " << spread_pairs.size() << " " << scores.size() << endl;
				cnt_it++;
				std::string IDD = createPairID(node1, sp1, node2, sp2);
				if(seedcheck.find(IDD) == seedcheck.end())
				{
					seedcheck.insert(IDD);
					for(int i=0; i < g[sp1]._nEdges[node1]; i++)
					{
						for(int j=0; j < g[sp2]._nEdges[node2]; j++)
						{
							int tnode1 = g[sp1]._edges[node1][i];
							int tnode2 = g[sp2]._edges[node2][j];							
							std::string ID1 = createID(tnode1, sp1);
							std::string ID2 = createID(tnode2, sp2);
							if(check_flag ||  (matched.find(ID1) == matched.end() && matched.find(ID2) == matched.end()))
							{		
								std::string ID = createPairID(tnode1, sp1, tnode2, sp2);
								if(score_map.find(ID) == score_map.end())
								{						
									score_map[ID] = 1;
									ProtPair pp; pp.prot1 = tnode1; pp.sp1 = sp1; pp.prot2 = tnode2; pp.sp2 = sp2; pp.value = 1;
									scores.insert(pp);
								}
								else
								{							
									int val = score_map[ID];
									ProtPair pp; pp.prot1 = tnode1; pp.sp1 = sp1; pp.prot2 = tnode2; pp.sp2 = sp2; pp.value = val;
									std::multiset<ProtPair, CompareProtPairs>::iterator itt = scores.find(pp);
									if(itt != scores.end())
										scores.erase(itt);
									val+=1;
									score_map[ID] = val;
									pp.value = val;
									scores.insert(pp);
								}	
							}
						}
					}
				}
			}	
		}*/
	}
	cout << "final score size: " << scores.size() << endl;
	cout << "size: " << tuples.size() << endl;
	return tuples;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::map <int, vector<int>> Alignment::triple(Graph g[], std::list<Triplet> & seed)
{
	std::map< std::string, int> score_map;
	std::multiset<ProtPair, CompareProtPairs> scores;
	std::set<std::string> seedcheck;
	int nb_networks = 3;
	std::set<int> * matched_nodes = new std::set<int> [nb_networks];//set of matched nodes
	int check = 0;
	int t_cnt = 0;
	int cnt = 0;
	int tuple[nb_networks];
	int nb_tuple = 0;
	std::map <int , vector<int> > tuples;
	
	std::map<int, int> * node_tuple = new std::map<int, int> [nb_networks];

	for ( list< Triplet>::iterator it = seed.begin(); it != seed.end(); ++it )
	{	
		Triplet q; q.node1 = it->node1; q.node2 = it->node2; q.node3 = it->node3; q.value = 3;
		for(int i = 0; i < nb_networks; i++)
			tuple[i] = -1;
		tuple[0] = it->node1; tuple[1] = it->node2; tuple[2] = it->node3;
		matched_nodes[0].insert(it->node1);
		matched_nodes[1].insert(it->node2);
		matched_nodes[2].insert(it->node3);
		std::vector<int> v;
		for(int i = 0; i < nb_networks; i++)
		{
			v.push_back(tuple[i]);
			if(tuple[i] != -1)
			{
				node_tuple[i][tuple[i]] = nb_tuple;
			}	
		}
		tuples[nb_tuple] = v;
		nb_tuple++;
		cnt += 1;
	}
	//return tuples;
	//////////////////////////////////////////////////////////////////////////
	for (list<Triplet>::iterator it = seed.begin(); it != seed.end(); ++it )
	{
		int node1 = it->node1; int node2 = it->node2; int node3 = it->node3;
		int node_array[nb_networks]; node_array[0] = node1; node_array[1] = node2;  node_array[2] = node3;
		for(int sp1 = 0; sp1 < nb_networks; sp1++)
		{
			for(int sp2 = sp1+1; sp2 < nb_networks; sp2++)
			{				
				if(node_array[sp1] != -1 and node_array[sp2] != -1)
				{
					std::string IDD = createPairID(node_array[sp1], sp1, node_array[sp2], sp2);
					if(seedcheck.find(IDD) == seedcheck.end())
					{
						seedcheck.insert(IDD);						
						for(int i=0; i < g[sp1]._nEdges[node_array[sp1]]; i++)
						{
							for(int j=0; j < g[sp2]._nEdges[node_array[sp2]]; j++)
							{
								//if(matched_nodes[sp1].find(g[sp1]._edges[node_array[sp1]][i]) == matched_nodes[sp1].end() && matched_nodes[sp2].find(g[sp2]._edges[node_array[sp2]][j]) == matched_nodes[sp2].end())
								{
									int node1 = g[sp1]._edges[node_array[sp1]][i];
									int node2 = g[sp2]._edges[node_array[sp2]][j];
									std::string ID = createPairID(node1, sp1, node2, sp2);
									if(score_map.find(ID) == score_map.end())
									{
										score_map[ID] = 1;
										ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 1;
										scores.insert(pp);
									}
									else
									{
										int val = score_map[ID];
										ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = val;
										std::multiset<ProtPair, CompareProtPairs>::iterator it = scores.find(pp);
										if(it != scores.end())
											scores.erase(it);
										val+=1;
										score_map[ID] = val;
										pp.value = val;
										scores.insert(pp);
									}	
								}
							}
						}
					}
				}
			}

		}
	}
	cout << "define2" << endl;
	/////////////////////////////////////////////////////////////////////////
	
	std::multiset<ProtPair, CompareProtPairs>::iterator it = scores.begin();
	cout << "score size: " << scores.size() << endl;
	bool expand = true;
	for(int run = 0; run < 1; run++)
	{
		while(scores.size() > 0 and it->value > 1)
		{	
		
			int node1 = it->prot1;
			int node2 = it->prot2;
			int sp1 = it->sp1;
			int sp2 = it->sp2;
		
			bool spread_flag = false;
			std::vector<ProtPair> spread_pairs;
			if(node_tuple[sp1].find(node1) == node_tuple[sp1].end() and node_tuple[sp2].find(node2) == node_tuple[sp2].end())
			{
				//if(node1 == node2)
					//cout << seedcheck.size() << " ++++ " << node1 << "\t" << node2 << "\t" << sp1 << "\t" << sp2 << endl;;			
				std::vector<int> v;
				for(int i = 0; i < nb_networks; i++)
				{
					v.push_back(-1);
				}
				if(node1 == -1 or node2 == -1)
				{
					cout << "hey" << endl;
				}
				v[sp1] = node1;
				v[sp2] = node2;
				node_tuple[sp1][node1] = nb_tuple;
				node_tuple[sp2][node2] = nb_tuple;
				tuples[nb_tuple] = v;
				matched_nodes[sp1].insert(node1);
				matched_nodes[sp2].insert(node2);
				nb_tuple++;
				spread_flag = true;
				ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 1;
				spread_pairs.push_back(pp);
			}

			if(node_tuple[sp1].find(node1) != node_tuple[sp1].end() and node_tuple[sp2].find(node2) == node_tuple[sp2].end())
			{
				int tuple_ind = node_tuple[sp1][node1];
				if(tuples[tuple_ind][sp2] == -1)
				{
					std::vector<int> v;
					v = tuples[tuple_ind];
					v[sp2] = node2;
					tuples[tuple_ind] = v;
					node_tuple[sp2][node2] = tuple_ind;
					matched_nodes[sp2].insert(node2);
					spread_flag = true;
					for(int i = 0; i < nb_networks; i++)
					{
						if(i != sp2 and tuples[tuple_ind][i] != -1)
						{
							ProtPair pp; pp.prot1 = tuples[tuple_ind][i]; pp.sp1 = i; pp.prot2 = node2; pp.sp2 = sp2; pp.value = 0;
							spread_pairs.push_back(pp);
						}
					}
				}
			}
		
			if(node_tuple[sp1].find(node1) == node_tuple[sp1].end() and node_tuple[sp2].find(node2) != node_tuple[sp2].end())
			{
				int tuple_ind = node_tuple[sp2][node2];
				if(tuples[tuple_ind][sp1] == -1)
				{
					std::vector<int> v;
					v = tuples[tuple_ind];
					v[sp1] = node1;
					tuples[tuple_ind] = v;
					node_tuple[sp1][node1] = tuple_ind;
					matched_nodes[sp1].insert(node1);
					spread_flag = true;
					for(int i = 0; i < nb_networks; i++)
					{
						if(i != sp1 and tuples[tuple_ind][i] != -1)
						{
							ProtPair pp; pp.prot1 = node1; pp.sp1 = sp1; pp.prot2 = tuples[tuple_ind][i]; pp.sp2 = i; pp.value = 0;
							spread_pairs.push_back(pp);
						}
					}
				}
			}
	
			//it = scores.erase(it);
			//std::string ID = createPairID(node1, sp1, node2, sp2);
			//score_map.erase(ID);
			/*if(spread_flag)
			{		
				for(std::vector<ProtPair>::iterator itv = spread_pairs.begin(); itv != spread_pairs.end(); ++itv)
				{
					int node1 = itv->prot1;
					int node2 = itv->prot2;
					int sp1 = itv->sp1;
					int sp2 = itv->sp2;
					for(int i=0; i < g[sp1]._nEdges[node1]; i++)
					{
						for(int j=0; j < g[sp2]._nEdges[node2]; j++)
						{
							if(matched_nodes[sp1].find(g[sp1]._edges[node1][i]) == matched_nodes[sp1].end() && matched_nodes[sp2].find(g[sp2]._edges[node2][j]) == matched_nodes[sp2].end())
							{
					
								int tnode1 = g[sp1]._edges[node1][i];
								int tnode2 = g[sp2]._edges[node2][j];
								std::string ID = createPairID(tnode1, sp1, tnode2, sp2);
								if(score_map.find(ID) == score_map.end())
								{						
									score_map[ID] = 1;
									ProtPair pp; pp.prot1 = tnode1; pp.sp1 = sp1; pp.prot2 = tnode2; pp.sp2 = sp2; pp.value = 1;
									scores.insert(pp);
								}
								else
								{							
									int val = score_map[ID];
									ProtPair pp; pp.prot1 = tnode1; pp.sp1 = sp1; pp.prot2 = tnode2; pp.sp2 = sp2; pp.value = val;
									std::multiset<ProtPair, CompareProtPairs>::iterator itt = scores.find(pp);
									if(itt != scores.end())
										scores.erase(itt);
									val+=1;
									score_map[ID] = val;
									pp.value = val;
									scores.insert(pp);
								}	
							}
						}
					}
				}
				it = scores.begin();
			}
			else
			{
				it = scores.erase(it);
				std::string ID = createPairID(node1, sp1, node2, sp2);
				score_map.erase(ID);
			}*/
		
			if(spread_flag)
			{		
				for(std::vector<ProtPair>::iterator itv = spread_pairs.begin(); itv != spread_pairs.end(); ++itv)
				{			
					int node1 = itv->prot1;
					int node2 = itv->prot2;
					int sp1 = itv->sp1;
					int sp2 = itv->sp2;
					std::string IDD = createPairID(node1, sp1, node2, sp2);
					if(seedcheck.find(IDD) == seedcheck.end())
					{
						seedcheck.insert(IDD);
						for(int i=0; i < g[sp1]._nEdges[node1]; i++)
						{
							for(int j=0; j < g[sp2]._nEdges[node2]; j++)
							{
								//if(matched_nodes[sp1].find(g[sp1]._edges[node1][i]) == matched_nodes[sp1].end() && matched_nodes[sp2].find(g[sp2]._edges[node2][j]) == matched_nodes[sp2].end())
								{
					
									int tnode1 = g[sp1]._edges[node1][i];
									int tnode2 = g[sp2]._edges[node2][j];
									std::string ID = createPairID(tnode1, sp1, tnode2, sp2);
									if(score_map.find(ID) == score_map.end())
									{						
										score_map[ID] = 1;
										ProtPair pp; pp.prot1 = tnode1; pp.sp1 = sp1; pp.prot2 = tnode2; pp.sp2 = sp2; pp.value = 1;
										scores.insert(pp);
									}
									else
									{							
										int val = score_map[ID];
										ProtPair pp; pp.prot1 = tnode1; pp.sp1 = sp1; pp.prot2 = tnode2; pp.sp2 = sp2; pp.value = val;
										std::multiset<ProtPair, CompareProtPairs>::iterator itt = scores.find(pp);
										if(itt != scores.end())
											scores.erase(itt);
										val+=1;
										score_map[ID] = val;
										pp.value = val;
										scores.insert(pp);
									}	
								}
							}
						}
					}
				}
				it = scores.begin();
			}
			else
			{
				it = scores.erase(it);
				std::string ID = createPairID(node1, sp1, node2, sp2);
				score_map.erase(ID);
			}
		}

	}
	cout << "final score size: " << scores.size() << endl;
	cout << "size: " << tuples.size() << endl;
	return tuples;
}

