#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <graph.h>
using namespace std;
#include "alignment.h"

using namespace std;

int main( int argc, char * argv[] )
{

	int nb_networks = atoi(argv[1]);
	string network_directory = argv[2];
	std::string sps[nb_networks];
	std::string graph_files[nb_networks];
	std::set<int> nodes [nb_networks];
	for(int net = 0; net < nb_networks; net++)
	{
		sps[net] = argv[3 + net];
		graph_files[net] = network_directory + "/" + argv[3 + net] + ".interaction";
		cout << sps[net] << endl;
		cout << graph_files[net] << endl;
	}

	std::string blastfile = argv[3 + nb_networks];
	std::string filename = argv[3 + nb_networks + 1];
	cout << blastfile << endl;
	cout << filename << endl;

	std::map<int , std::string> * inttonodeID = new std::map<int , std::string> [nb_networks];
	std::map<std::string , int> * nodeIDtoint = new std::map<std::string , int> [nb_networks];
	Graph * g_array = new Graph [nb_networks];
	for(int i = 0; i < nb_networks; i++)
	{
		Graph * sg = new Graph(false);
		sg->readGraph( graph_files[i]);
		g_array[i] = *sg;
		cout << g_array[i].getNNodes() << endl;
		inttonodeID[i] = g_array[i]._inttonodeID;
		nodeIDtoint[i] = g_array[i]._nodeIDtoint;
	}

	for(int i = 0; i < nb_networks; i++){

		cout << "Number of nodes: " << g_array[i].getNNodes() << endl;
	}



	std::list< Match>** seeds_pairs = new std::list< Match>*[nb_networks];
	for (int i = 0; i < nb_networks; ++i)
	{
	   seeds_pairs[i] = new std::list< Match> [nb_networks]; //pairwise alignments
	}
	std::map <int , vector<int> > seeds;
	std::set<std::string> seedcheck;
	Match m;
	cout << "Seeds!" << endl;
	int seed_ind = 0;


	std::ifstream infile(blastfile);

	std::string line;
	while (std::getline(infile, line))
	{
	    std::istringstream iss(line);
	    vector<string> tokens{istream_iterator<string>{iss},
                      istream_iterator<string>{}};
        if (tokens.size() != nb_networks)
        {

        	cout << "The seedtuple file is erroneous!" << endl;
        	return 0;
        }
        std::vector<int> v;
		for(int i = 0; i < nb_networks; i++)
		{
			if (tokens[i] != "-1")
			{
				v.push_back(nodeIDtoint[i][tokens[i]]);
			}
			else
			{
				v.push_back(-1);
			}
		}
		seeds[seed_ind] = v;
		seed_ind++;		
	}

	
	

	Alignment * align = new Alignment();		

	std::map <int , vector<int> > tuples = align->multipleK(g_array, nb_networks, seeds, true);
	cout << "end!" << endl;
	ofstream output( filename.c_str() );
	cout << "Alignment size: " << tuples.size() << endl; 
	for(int i = 0; i < tuples.size(); i++)
	{
		std::string match_line = "";
		int t_node = 0;
		for(int j = 0; j < nb_networks; j++)
		{
			if(tuples[i][j] == -1)
			{
				match_line += "-1";
			}
			else
			{
				match_line += inttonodeID[j][tuples[i][j]];
				t_node += 1;
			}
			if(j != nb_networks - 1)
			{
				match_line+="\t";
			}
		}
		if(t_node > 1)
		{ 
			output << match_line << std::endl;
		}

	}
	output.close();
	
	return 0;
}


