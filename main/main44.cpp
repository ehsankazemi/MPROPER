#include "rg.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <time.h>
#include <algorithm>
enum Type { ErdosRenyi, BarabasiAlbert, WattsStrogatz }; //0 for GNP, 1 for BA and 2 for WS

using namespace std;

std::string createMatchID( int i, int j )
{
    std::stringstream sstm;
    sstm << i << "$" << j;
    return sstm.str();
}

void untokenizeMatchID( std::string str, int & first, int & second )
{
    std::string ES = EdgeSeparator;
    int index = str.find( ES );
    first  = atoi(str.substr( 0, index ).c_str());
    if ( index != -1 )
        second = atoi(str.substr( str.find( ES ) + ES.size() ).c_str());
}

std::vector < std::map< std::string, int> > blast_scores(std::string sp1, std::string sp2, int min_blast)
{
	std::vector < std::map< std::string, int> > blasts;
	std::map< std::string, int> lrblast;
	std::map< std::string, int> rlblast;
	int u, v;
	int n;
	string s1,s2,s3,s4;
    Match m;
	string blastfile = "proc/sorted_blast_"+ sp1 +"-" + sp2 +".blast";
	
	ifstream file;
	file.open( blastfile.c_str() );
	int cnt = 0;
	while (file >> s1 >> u >> s2 >> s3  >> v >> s4 >> n  )
	{		
		cnt++;
		if(n < min_blast)
			break;
		m.lnode = u; m.rnode = v; m.value = n;
		lrblast[createMatchID(u,v)] = n;
		rlblast[createMatchID(v,u)] = n;
	}
	blasts.push_back(lrblast);
	blasts.push_back(rlblast);
	cout << blastfile << "\t" << lrblast.size() << endl;
	return blasts;
}




std::vector <std::map< int, std::list<NodeVal> > > read_blast_scores(std::string sp1, std::string sp2, int min_blast)
{
	std::vector <std::map< int, std::list<NodeVal> > > blasts;
	std::map< int, std::list<NodeVal> > lrblast;
	std::map< int, std::list<NodeVal> > rlblast;
	int u, v;
	int n;
	string s1,s2,s3,s4;
    NodeVal nvlr;
	NodeVal nvrl;
	string blastfile = "proc/sorted_blast_"+ sp1 +"-" + sp2 +".blast";
	
	ifstream file;
	file.open( blastfile.c_str() );
	int cnt = 0;
	while (file >> s1 >> u >> s2 >> s3  >> v >> s4 >> n  )
	{		
		cnt++;
		if(n < min_blast)
			break;
		nvlr.rnode = v; nvlr.value = n;
		nvrl.rnode = u; nvrl.value = n;
		lrblast[u].push_back(nvlr);
		rlblast[v].push_back(nvrl);
	}
		
	cout << blastfile << "\t" << lrblast.size() << endl;
	blasts.push_back(lrblast);
	blasts.push_back(rlblast);
	return blasts;
}


std::vector <std::map< int, std::list<NodeVal> > > reread_blast_scores(std::string sp1, std::string sp2, int min_blast, std::set<int> nodes1, std::set<int> nodes2 )
{
	std::vector <std::map< int, std::list<NodeVal> > > blasts;
	std::map< int, std::list<NodeVal> > lrblast;
	std::map< int, std::list<NodeVal> > rlblast;
	int u, v;
	int n;
	string s1,s2,s3,s4;
    NodeVal nvlr;
	NodeVal nvrl;
	string blastfile = "proc/sorted_blast_"+ sp1 +"-" + sp2 +".blast";
	
	ifstream file;
	file.open( blastfile.c_str() );
	int cnt = 0;
	while (file >> s1 >> u >> s2 >> s3  >> v >> s4 >> n  )
	{		
		cnt++;
		if(n < min_blast)
			break;
		if(nodes1.find(u) == nodes1.end() and nodes2.find(v) == nodes2.end())
		{
			nvlr.rnode = v; nvlr.value = n;
			nvrl.rnode = u; nvrl.value = n;
			lrblast[u].push_back(nvlr);
			rlblast[v].push_back(nvrl);
		}
	}
		
	cout << blastfile << "\t" << lrblast.size() << endl;
	blasts.push_back(lrblast);
	blasts.push_back(rlblast);
	return blasts;
}

void load_blasts(std::map< int, std::list<NodeVal> >** blasts, std::string sps[5], int min_blast, std::set<int> nodes [5])
{
	for(int i = 0; i < 5; i++)
	{
		for(int j = i+1; j < 5; j++)
		{
			if(i != j)
			{
				std::vector <std::map< int, std::list<NodeVal> > > blast_vecs;
				blast_vecs = reread_blast_scores(sps[i], sps[j], min_blast, nodes[i], nodes[j]);
				blasts[i][j] = blast_vecs.at(0);
				blasts[j][i] = blast_vecs.at(1);
			}
		}
	}
}

int main( int argc, char * argv[] )
{
	srand( time( NULL ) );
	int min_blast = atoi( argv[1] );
	string sp1 = "ce";
	string graph_file1 = "proc/" + sp1 + ".txt";
	string sp2 = "dm";
	string graph_file2 = "proc/" + sp2 + ".txt";
	string sp3 = "hs";
	string graph_file3 = "proc/" + sp3 + ".txt";
	string sp4 = "mm";
	string graph_file4 = "proc/" + sp4 + ".txt";
	string sp5 = "sc";
	string graph_file5 = "proc/" + sp5 + ".txt";

  	
  	std::set<int> nodes [5];
	std::string sps[5];
	sps[0] = sp1;
	sps[1] = sp2;
	sps[2] = sp3;
	sps[3] = sp4;
	sps[4] = sp5;

	
string blastfile = "scripts/hh_"+ sp1 + "_" +sp2 + "_" + sp3 + "_" + sp4 + "_" + sp5  + "_" + std::string(argv[1]) + ".txt";
	ofstream output( blastfile.c_str() );
	bool flag;
	int cnt = 0;
	clock_t tStart = clock();
	
	int myints[] = {0, 1, 2, 3, 4};

 	int	sub_sps_ind [5] = {0, 1, 2, 3, 4};


	int	sps_ind [5] = {0, 1, 2, 3, 4};
	int thr;
	int up_bound = 200;
	int run = (up_bound - min_blast) / 50;
	int nb = 2;
	int dist = 100;
	for(int r = nb; r >= 0; r--)
	{
		thr = min_blast + dist * r;
		std::map< int, std::list<NodeVal> >** blasts =new std::map< int, std::list<NodeVal> >*[5];	
		for (int i = 0; i < 5; ++i)
		{
		   blasts[i] = new std::map< int, std::list<NodeVal> >[5];
		}
		load_blasts(blasts, sps, thr, nodes);
		do
		{
			
			for(std::map< int, std::list<NodeVal> >::iterator it1 = blasts[sps_ind[0]][sps_ind[1]].begin(); it1 != blasts[sps_ind[0]][sps_ind[1]].end(); it1++)
			{
				flag = false;
				int prot1 = it1->first;
				if(nodes[sps_ind[0]].find(prot1) == nodes[sps_ind[0]].end())
				{
					std::list<NodeVal>::iterator itt1 = blasts[sps_ind[0]][sps_ind[1]][prot1].begin();
					while(itt1 != blasts[sps_ind[0]][sps_ind[1]][prot1].end()){
						if(flag)
							break;
						int prot2 = itt1->rnode;
						int val12 = itt1->value;
						if(val12 >= thr)
						{
							if(nodes[sps_ind[1]].find(prot2) == nodes[sps_ind[1]].end()){
								if(blasts[sps_ind[1]][sps_ind[2]].find(prot2) != blasts[sps_ind[1]][sps_ind[2]].end()){
									std::list<NodeVal>::iterator itt2 = blasts[sps_ind[1]][sps_ind[2]][prot2].begin();
									while(itt2 != blasts[sps_ind[1]][sps_ind[2]][prot2].end()){
										if(flag) break;						
										int prot3 = itt2->rnode;
										int val23 = itt2->value;
										if(val23 >= thr)
										{
											if(nodes[sps_ind[2]].find(prot3) == nodes[sps_ind[2]].end()){
												if(blasts[sps_ind[2]][sps_ind[3]].find(prot3) != blasts[sps_ind[2]][sps_ind[3]].end()){
													std::list<NodeVal>::iterator itt3 = blasts[sps_ind[2]][sps_ind[3]][prot3].begin();
													while(itt3 != blasts[sps_ind[2]][sps_ind[3]][prot3].end()){
														if(flag) break;						
														int prot4 = itt3->rnode;
														int val34 = itt3->value;
														if(val34 >= thr)
														{
															if(nodes[sps_ind[3]].find(prot4) == nodes[sps_ind[3]].end()){
																if(blasts[sps_ind[3]][sps_ind[4]].find(prot4) != blasts[sps_ind[3]][sps_ind[4]].end()){
																	std::list<NodeVal>::iterator itt4 = blasts[sps_ind[3]][sps_ind[4]][prot4].begin();
																	while(itt4 != blasts[sps_ind[3]][sps_ind[4]][prot4].end()){
																		int prot5 = itt4->rnode;
																		int val45 = itt4->value;
																		if(val45 >= thr)
																		{				
																			if(nodes[sps_ind[4]].find(prot5) == nodes[sps_ind[4]].end()){
																				cnt++;
																				int prots [5] = {-1, -1, -1, -1, -1};
																				prots[sps_ind[0]] = prot1;
																				prots[sps_ind[1]] = prot2;
																				prots[sps_ind[2]] = prot3;
																				prots[sps_ind[3]] = prot4;
																				prots[sps_ind[4]] = prot5;

																				nodes[sps_ind[0]].insert(prot1);
																				nodes[sps_ind[1]].insert(prot2);
																				nodes[sps_ind[2]].insert(prot3);
																				nodes[sps_ind[3]].insert(prot4);
																				nodes[sps_ind[4]].insert(prot5);
																				output << prots[0] << "\t" << prots[1] << "\t" << prots[2] << "\t" << prots[3] << "\t" << prots[4] << std::endl;
																				cout << cnt << endl;
																				std::map< int, std::list<NodeVal> >::iterator it = blasts[sps_ind[1]][sps_ind[2]].find(prot2);
																				if(it != blasts[sps_ind[1]][sps_ind[2]].end()) blasts[sps_ind[1]][sps_ind[2]].erase(it);

																				it = blasts[sps_ind[2]][sps_ind[3]].find(prot3); if(it != blasts[sps_ind[2]][sps_ind[3]].end()) blasts[sps_ind[2]][sps_ind[3]].erase(it);

																				it = blasts[sps_ind[3]][sps_ind[4]].find(prot4);
																				if(it != blasts[sps_ind[3]][sps_ind[4]].end()) blasts[sps_ind[3]][sps_ind[4]].erase(it);
																				flag = true;
																				break;
																			}
																			itt4 = blasts[sps_ind[3]][sps_ind[4]][prot4].erase(itt4);
																		}
																		else itt4++;
																	}
																}
																itt3++;
															}
															else itt3 = blasts[sps_ind[2]][sps_ind[3]][prot3].erase(itt3);
														}
														else itt3++;
													}
												}
												itt2++;
											}
											else itt2 = blasts[sps_ind[1]][sps_ind[2]][prot2].erase(itt2);
										}
										else itt2++;
									}
								}
								itt1++;
							}
							else itt1 = blasts[sps_ind[0]][sps_ind[1]][prot1].erase(itt1);
						}
						else itt1++;
					}
				}
			//cout << "+++++++++++++++++" << endl;
			}
		} while ( std::next_permutation(sps_ind, sps_ind + 5) );
	for(int i = 0; i < 5; ++i) {
		delete [] blasts[i];
	}
	delete [] blasts;
	}
	cout << "Done!" << endl;
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(int r = nb; r >= 0; r--)
	{
		thr = min_blast + dist * r;
		std::map< int, std::list<NodeVal> >** blasts =new std::map< int, std::list<NodeVal> >*[5];	
		for (int i = 0; i < 5; ++i)
		{
		   blasts[i] = new std::map< int, std::list<NodeVal> >[5];
		}
		load_blasts(blasts, sps, thr, nodes);
		do
		{
			cout << "+++++++++++++++++++++" << endl;
			for(std::map< int, std::list<NodeVal> >::iterator it1 = blasts[sub_sps_ind[0]][sub_sps_ind[1]].begin(); it1 != blasts[sub_sps_ind[0]][sub_sps_ind[1]].end(); it1++)
			{
				flag = false;
				int prot1 = it1->first;
				if(nodes[sub_sps_ind[0]].find(prot1) == nodes[sub_sps_ind[0]].end())
				{
					std::list<NodeVal>::iterator itt1 = blasts[sub_sps_ind[0]][sub_sps_ind[1]][prot1].begin();
					while(itt1 != blasts[sub_sps_ind[0]][sub_sps_ind[1]][prot1].end()){
						if(flag)
							break;
						int prot2 = itt1->rnode;
						int val12 = itt1->value;	
						if(nodes[sub_sps_ind[1]].find(prot2) == nodes[sub_sps_ind[1]].end()){
							if(blasts[sub_sps_ind[1]][sub_sps_ind[2]].find(prot2) != blasts[sub_sps_ind[1]][sub_sps_ind[2]].end()){
								std::list<NodeVal>::iterator itt2 = blasts[sub_sps_ind[1]][sub_sps_ind[2]][prot2].begin();
								while(itt2 != blasts[sub_sps_ind[1]][sub_sps_ind[2]][prot2].end()){
									if(flag) break;						
									int prot3 = itt2->rnode;
									int val23 = itt2->value;
									if(nodes[sub_sps_ind[2]].find(prot3) == nodes[sub_sps_ind[2]].end()){
										if(blasts[sub_sps_ind[2]][sub_sps_ind[3]].find(prot3) != blasts[sub_sps_ind[2]][sub_sps_ind[3]].end()){
											std::list<NodeVal>::iterator itt3 = blasts[sub_sps_ind[2]][sub_sps_ind[3]][prot3].begin();
											while(itt3 != blasts[sub_sps_ind[2]][sub_sps_ind[3]][prot3].end()){
												int prot4 = itt3->rnode;
												int val34 = itt3->value;
												if(nodes[sub_sps_ind[3]].find(prot3) == nodes[sub_sps_ind[3]].end()){
													cnt++;
													int prots [5] = {-1, -1, -1, -1, -1};
													prots[sub_sps_ind[0]] = prot1;
													prots[sub_sps_ind[1]] = prot2;
													prots[sub_sps_ind[2]] = prot3;
													prots[sub_sps_ind[3]] = prot4;

													nodes[sub_sps_ind[0]].insert(prot1);
													nodes[sub_sps_ind[1]].insert(prot2);
													nodes[sub_sps_ind[2]].insert(prot3);
													nodes[sub_sps_ind[3]].insert(prot4);

													output << prots[0] << "\t" << prots[1] << "\t" << prots[2] << "\t" << prots[3] << "\t" << prots[4] << std::endl;
													cout << cnt << endl;
													std::map< int, std::list<NodeVal> >::iterator it = blasts[sub_sps_ind[1]][sub_sps_ind[2]].find(prot2);
													if(it != blasts[sub_sps_ind[1]][sub_sps_ind[2]].end()) blasts[sub_sps_ind[1]][sub_sps_ind[2]].erase(it);

													it = blasts[sub_sps_ind[2]][sub_sps_ind[3]].find(prot3); if(it != blasts[sub_sps_ind[2]][sub_sps_ind[3]].end()) blasts[sub_sps_ind[2]][sub_sps_ind[3]].erase(it);
													flag = true;
													break;
												}
												itt3 = blasts[sub_sps_ind[2]][sub_sps_ind[3]][prot3].erase(itt3);
											}
										}
										itt2++;
									}
									else itt2 = blasts[sub_sps_ind[1]][sub_sps_ind[2]][prot2].erase(itt2);
								}
							}
							itt1++;
						}
						else itt1 = blasts[sub_sps_ind[0]][sub_sps_ind[1]][prot1].erase(itt1);
					}
				}
			}
		} while ( std::next_permutation(sub_sps_ind, sub_sps_ind + 5) );
	for(int i = 0; i < 5; ++i) {
		delete [] blasts[i];
	}
	delete [] blasts;		
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	//min_blast = 40;
	sub_sps_ind[0] = 0;
	sub_sps_ind[1] = 1;
	sub_sps_ind[2] = 2;
	sub_sps_ind[3] = 3;
	sub_sps_ind[4] = 4;
	for(int r = nb; r >= 0; r--)
	{
		thr = min_blast + dist * r;
		std::map< int, std::list<NodeVal> >** blasts =new std::map< int, std::list<NodeVal> >*[5];	
		for (int i = 0; i < 5; ++i)
		{
		   blasts[i] = new std::map< int, std::list<NodeVal> >[5];
		}
		load_blasts(blasts, sps, thr, nodes);
		do
		{
			cout << "+++++++++++++++++++++" << endl;
			for(std::map< int, std::list<NodeVal> >::iterator it1 = blasts[sub_sps_ind[0]][sub_sps_ind[1]].begin(); it1 != blasts[sub_sps_ind[0]][sub_sps_ind[1]].end(); it1++)
			{
				flag = false;
				int prot1 = it1->first;
				if(nodes[sub_sps_ind[0]].find(prot1) == nodes[sub_sps_ind[0]].end())
				{
					std::list<NodeVal>::iterator itt1 = blasts[sub_sps_ind[0]][sub_sps_ind[1]][prot1].begin();
					while(itt1 != blasts[sub_sps_ind[0]][sub_sps_ind[1]][prot1].end()){
						if(flag)
							break;
						int prot2 = itt1->rnode;
						int val12 = itt1->value;	
						if(nodes[sub_sps_ind[1]].find(prot2) == nodes[sub_sps_ind[1]].end()){
							if(blasts[sub_sps_ind[1]][sub_sps_ind[2]].find(prot2) != blasts[sub_sps_ind[1]][sub_sps_ind[2]].end()){
								std::list<NodeVal>::iterator itt2 = blasts[sub_sps_ind[1]][sub_sps_ind[2]][prot2].begin();
								while(itt2 != blasts[sub_sps_ind[1]][sub_sps_ind[2]][prot2].end()){
									int prot3 = itt2->rnode;
									int val23 = itt2->value;
									if(nodes[sub_sps_ind[2]].find(prot3) == nodes[sub_sps_ind[2]].end()){
										cnt++;
										int prots [5] = {-1, -1, -1, -1, -1};
										prots[sub_sps_ind[0]] = prot1;
										prots[sub_sps_ind[1]] = prot2;
										prots[sub_sps_ind[2]] = prot3;


										nodes[sub_sps_ind[0]].insert(prot1);
										nodes[sub_sps_ind[1]].insert(prot2);
										nodes[sub_sps_ind[2]].insert(prot3);

										output << prots[0] << "\t" << prots[1] << "\t" << prots[2] << "\t" << prots[3] << "\t" << prots[4] << std::endl;
										cout << cnt << endl;
										std::map< int, std::list<NodeVal> >::iterator it = blasts[sub_sps_ind[1]][sub_sps_ind[2]].find(prot2);
										if(it != blasts[sub_sps_ind[1]][sub_sps_ind[2]].end()) blasts[sub_sps_ind[1]][sub_sps_ind[2]].erase(it);

										flag = true;
										break;
									}
									itt2 = blasts[sub_sps_ind[1]][sub_sps_ind[2]][prot2].erase(itt2);
								}
							}
							itt1++;
						}
						else itt1 = blasts[sub_sps_ind[0]][sub_sps_ind[1]][prot1].erase(itt1);
					}
				}
			}
		} while ( std::next_permutation(sub_sps_ind, sub_sps_ind + 5) );
	for(int i = 0; i < 5; ++i) {
		delete [] blasts[i];
	}
	delete [] blasts;	
	}
	//}
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	sub_sps_ind[0] = 0;
	sub_sps_ind[1] = 1;
	sub_sps_ind[2] = 2;
	sub_sps_ind[3] = 3;
	sub_sps_ind[4] = 4;
	std::map< int, std::list<NodeVal> >** blasts =new std::map< int, std::list<NodeVal> >*[5];	
	for (int i = 0; i < 5; ++i)
	{
	   blasts[i] = new std::map< int, std::list<NodeVal> >[5];
	}
	load_blasts(blasts, sps, min_blast, nodes);
	do
	{
		cout << "+++++++++++++++++++++" << endl;
		for(std::map< int, std::list<NodeVal> >::iterator it1 = blasts[sub_sps_ind[0]][sub_sps_ind[1]].begin(); it1 != blasts[sub_sps_ind[0]][sub_sps_ind[1]].end(); it1++)
		{
			flag = false;
			int prot1 = it1->first;
			if(nodes[sub_sps_ind[0]].find(prot1) == nodes[sub_sps_ind[0]].end())
			{
				std::list<NodeVal>::iterator itt1 = blasts[sub_sps_ind[0]][sub_sps_ind[1]][prot1].begin();
				while(itt1 != blasts[sub_sps_ind[0]][sub_sps_ind[1]][prot1].end()){
					int prot2 = itt1->rnode;
					int val12 = itt1->value;	
					if(nodes[sub_sps_ind[1]].find(prot2) == nodes[sub_sps_ind[1]].end()){
						cnt++;
						int prots [5] = {-1, -1, -1, -1, -1};
						prots[sub_sps_ind[0]] = prot1;
						prots[sub_sps_ind[1]] = prot2;



						nodes[sub_sps_ind[0]].insert(prot1);
						nodes[sub_sps_ind[1]].insert(prot2);


						output << prots[0] << "\t" << prots[1] << "\t" << prots[2] << "\t" << prots[3] << "\t" << prots[4] << std::endl;
						cout << cnt << endl;
						flag = true;
						break;
					}
					itt1 = blasts[sub_sps_ind[0]][sub_sps_ind[1]][prot1].erase(itt1);
				}
			}
		}
	} while ( std::next_permutation(sub_sps_ind, sub_sps_ind + 5) );
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////

	
	printf("Time taken: %.4fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
	output.close();
	for(int i = 0; i < 5; ++i) {
		delete [] blasts[i];
	}
	delete [] blasts;	
	return 0;
}
