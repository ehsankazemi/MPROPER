#include "rg.h"
#include "alignment.h"
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <ctime>

enum Type { ErdosRenyi, BarabasiAlbert, WattsStrogatz }; //0 for GNP, 1 for BA and 2 for WS

using namespace std;


int main( int argc, char * argv[] )
{
	srand( time( NULL ) );
	int N = atoi(argv[1]);
	int nb_networks = atoi(argv[2]);
	int nb_seeds = atoi(argv[3]);
	double ps = atof(argv[4]);
	double pt = atof(argv[5]);
	string filename = argv[ 6 ];
	ofstream output( filename.c_str() );

	int avg = 20;
	int k = 2;
	Graph * g;
	g = new Erdos_Renyi( N );
	reinterpret_cast< Erdos_Renyi * >( g )->build( avg, k ); //avg degree value and k-size

	Graph * g_array = new Graph [nb_networks];
	std::map<int , int> * inttonodeID = new std::map<int , int> [nb_networks];
	std::map<int , int> * nodeIDtoint = new std::map<int , int> [nb_networks];
	for(int i = 0; i < nb_networks; i++)
	{
		Graph * sg = new Graph(*g, ps, pt);
		g_array[i] = *sg;
		inttonodeID[i] = g_array[i]._inttonodeID;
		nodeIDtoint[i] = g_array[i]._nodeIDtoint;
	}
	delete g;


	std::map <int , vector<int> >** seeds_pairs = new std::map <int , vector<int> >*[nb_networks];
	for (int i = 0; i < nb_networks; ++i)
	{
	   seeds_pairs[i] = new std::map <int , vector<int> >[nb_networks]; //pairwise alignments
	}
	std::map <int , vector<int> > seeds;
	std::set<std::string> seedcheck;
	Match m;
	cout << "Seeds!" << endl;
	int seed_ind = 0;
	for(int i=0; i < nb_seeds; i++)
	{
		bool flag = true;
		
		int check = 0;
		while(flag)
		{
			int node = rand() % N;
			//cout << "Seeds0!" << endl;
			std::stringstream sstm;
			sstm << node << "$" << node;
			std::string ID = sstm.str();
			bool seed_flag = true;
			if(seedcheck.find(ID) != seedcheck.end())
				seed_flag = false;
			for(int j = 0; j < nb_networks; j++)
			{
				if(nodeIDtoint[j].find(node) == nodeIDtoint[j].end())
					seed_flag = false;
			}
			if(seed_flag)
			{
				std::vector<int> v;
				for(int i = 0; i < nb_networks; i++)
				{
					v.push_back(nodeIDtoint[i][node]);
				}		
				seeds[seed_ind] = v;
				
				seedcheck.insert(ID);
				flag = false;
				for(int ind1 = 0; ind1 < nb_networks; ind1++)
				{
					for(int ind2 = ind1+1; ind2 < nb_networks; ind2++)
					{
						std::vector<int> v;
						v.push_back(nodeIDtoint[ind1][node]);
						v.push_back(nodeIDtoint[ind2][node]);
						m.lnode = nodeIDtoint[ind1][node]; m.rnode = nodeIDtoint[ind2][node]; m.value = 1;
						seeds_pairs[ind1][ind2][seed_ind] = v;
					}
				}
				seed_ind++;
			}
		}
	}
	
	

	Alignment * align = new Alignment();		
	Graph * g_p = new Graph [2];		
	for(int ind1 = 0; ind1 < nb_networks; ind1++)
	{
		for(int ind2 = ind1+1; ind2 < nb_networks; ind2++)
		{
			g_p[0] = g_array[ind1];
			g_p[1] = g_array[ind2];
			std::map <int , vector<int> > tuples = align->multipleK(g_p, 2, seeds_pairs[ind1][ind2], true);
			int c = 0;
			int w = 0;
			for(int i = 0; i < tuples.size(); i++)
			{
				if(inttonodeID[ind1][tuples[i][0]] == inttonodeID[ind2][tuples[i][1]])
					c++;
				else
					w++;
			}
			cout << ind1 << " - " << ind2 << " +++ " << c << "\t" << w << endl;
			output << ind1 << "\t" << ind2 << "\t" << c << "\t" << w << std::endl;
		}
	}
	std::map <int , vector<int> > tuples = align->multipleK(g_array, nb_networks, seeds, true);
	cout << "end!" << endl;
	int c = 0;
	int w = 0;
	int gc = 0;

	int count[nb_networks + 1];
	for(int i = 0; i < nb_networks + 1; i++)
	{
		count[i] = 0;
	}

	for(int i = 0; i < tuples.size(); i++)
	{
		bool c_flag = true;
		int sp_cnt = 0;
		int t_node = -1;
		for(int j = 0; j < nb_networks; j++)
		{
			if(t_node == -1)
			{
				if(tuples[i][j] != -1)
				{
					t_node = inttonodeID[j][tuples[i][j]];
					sp_cnt++;
				}
			}
			else
			{
				if(tuples[i][j] != -1)
				{
					if(inttonodeID[j][tuples[i][j]] == t_node)
					{
						sp_cnt++;
					}
					else
					{
						c_flag = false;
					}
				}
			}
		}
		if(t_node != -1)
		{
			if(c_flag)
			{
				count[sp_cnt]++;
			}
			else
			{
				w++;
			}
		}
	}
	cout << "Tuple size: " << tuples.size() << endl; 
	cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
	output << "multiple" << std::endl;
	for(int i = nb_networks; i > 1; i--)
	{
		cout << i << ": " << count[i] << endl;
		output << i << "\t" << count[i] << std::endl;
	}
	cout << "Wrong: "  << w << endl;
	output << w << std::endl;
	output.close();
	return 0;
}

