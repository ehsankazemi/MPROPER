from collections import defaultdict
import math
import pylab
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import community
import copy
import os
import sys
import time
import subprocess
import pickle
import random
from random import shuffle
def intersect(a, b):
     return list(set(a) & set(b))

def union(a, b):
     return list(set(a) | set(b))

def load_data(sp1, sp2):
	filename = '/home/90days/ekazemi/blast/proc/'
	dict1_prot_c = pickle.load(open(filename + sp1+'_prot_to_c.dat', 'rb'))
	dict1_c_prot = pickle.load(open(filename + sp1+'_c_to_prot.dat', 'rb'))
	dict2_prot_c = pickle.load(open(filename + sp2+'_prot_to_c.dat', 'rb'))
	dict2_c_prot = pickle.load(open(filename + sp2+'_c_to_prot.dat', 'rb'))
	
	f = open( filename + '../filtered/'+sp1+'.tab', 'r' )
	G1 = nx.Graph()
	cnt = 0
	for line in f:
		if line.find( '#' ) != -1: continue
		cnt+=1
		if cnt == 1: continue
		parsed = line.replace( "\n", "" ).split()
		if parsed[0] != parsed[1]:
			G1.add_edge(parsed[0],parsed[1])
	f = open( filename + '../filtered/'+sp2+'.tab', 'r' )
	G2 = nx.Graph()
	cnt = 0
	for line in f:
		if line.find( '#' ) != -1: continue
		cnt+=1
		if cnt == 1: continue
		parsed = line.replace( "\n", "" ).split()
		if parsed[0] != parsed[1]:
			G2.add_edge(parsed[0],parsed[1])
	return G1, G2, dict1_prot_c, dict1_c_prot, dict2_prot_c, dict2_c_prot


def load_graph(sp1):
	filename = '/home/90days/ekazemi/blast/proc/'
	dict1_prot_c = pickle.load(open(filename + sp1+'_prot_to_c.dat', 'rb'))
	dict1_c_prot = pickle.load(open(filename + sp1+'_c_to_prot.dat', 'rb'))
	
	f = open( filename + '../filtered/'+sp1+'.tab', 'r' )
	G1 = nx.Graph()
	cnt = 0
	for line in f:
		if line.find( '#' ) != -1: continue
		cnt+=1
		if cnt == 1: continue
		parsed = line.replace( "\n", "" ).split()
		if parsed[0] != parsed[1]:
			G1.add_edge(parsed[0],parsed[1])
	
	return G1, dict1_prot_c, dict1_c_prot

def community_detection(G):
	print nx.number_connected_components(G)
	G0 = nx.connected_component_subgraphs(G)[0]
	print len(G.nodes()), len(G0.nodes())
	partition = community.best_partition(G0)
	communities = defaultdict(list)
	for node in partition:
		communities[partition[node]].append(node)
	return partition, communities, G0

def dict_to_list(dict_list):
	list_list = []
	cnt = 0 
	for l in dict_list:
		list_list.append([])
		for e in dict_list[l]:
			list_list[cnt].append(e)
		cnt +=1
	return list_list

def goc_exp(matching):
	f = open( '../go_exp.txt', 'r' )
	dict_prot_go_all = defaultdict(list)

	for line in f:
		if line.find('GO') != -1:
		
			parsed = line.replace( "\n", "" ).split('|')
			for pars in parsed:
				ppars = pars.split()
				for p in ppars:
					if p.find('GO') != -1:
						dict_prot_go_all[parsed[0]].append(p)

	tot = 0
	go =0
	gocorrect = 0
	for node1 in matching:
		prot1 = node1
		prot2 = matching[node1]
		tot+=1
		if prot1 in dict_prot_go_all and prot2 in dict_prot_go_all:
			if len(dict_prot_go_all[prot1]) > 0 and len(dict_prot_go_all[prot2]) > 0:
				go+=1
				intsec = intersect(dict_prot_go_all[prot1], dict_prot_go_all[prot2])
				un = union(dict_prot_go_all[prot1], dict_prot_go_all[prot2])
				if len(intsec) > 0:
					gocorrect+=(len(intsec) * 1.0 / len(un))
	#gocorrect += 1

	#print tot
	#print go
	return gocorrect


def goc_all(sp1, sp2, matching, name):
	#GO1 = pickle.load(open(sp1+'_prot_'+ name + 'go.pkl', 'rb'))
	#GO2 = pickle.load(open(sp2+'_prot_'+ name + 'go.pkl', 'rb'))
	GO1 = pickle.load(open("GO"+ name +".dat", "rb"))
	GO2 = GO1
	filename = '/home/90days/ekazemi/blast/proc/'
	dict1_c_prot = pickle.load(open(filename + sp1+'_c_to_prot.dat', 'rb'))
	dict2_c_prot = pickle.load(open(filename + sp2+'_c_to_prot.dat', 'rb'))

	go =0
	gocorrect = 0
	tot = 0
	for pair in matching:
		prot1 = dict1_c_prot[int(pair)]
		prot2 = dict2_c_prot[int(matching[pair])]
		tot+= 1
		if prot1 in GO1 and prot2 in GO2:
			if len(GO1[prot1]) > 0 and len(GO2[prot2]) > 0:
				go+=1
				intsec = intersect(GO1[prot1], GO2[prot2])
				un = union(GO1[prot1], GO2[prot2])
				if len(intsec) > 0:
					gocorrect+=(len(intsec) * 1.0 / len(un))
	return gocorrect

def ICS(G1, G2, dict_match):
	size = len(dict_match)
	nodes = []
	nodes2 = []
	for node in dict_match:
		nodes.append(node)
		nodes2.append(dict_match[node])
	H = G1.subgraph(nodes)
	H2 = G2.subgraph(nodes2)
	tot_edge = 0
	correct_edge = 0
	G12 = nx.Graph()
	for node in nodes:
		G12.add_node(node)
	for e in H.edges():
		if e[0] in dict_match and e[1] in dict_match:
			tot_edge += 1
			if G2.has_edge(dict_match[e[0]],dict_match[e[1]]):
				G12.add_edge(e[0],e[1])
				correct_edge += 1
	tot_edge2 = 0
	for e in H2.edges():
		if e[0] in nodes2 and e[1] in nodes2:
			tot_edge2 += 1

	tot = tot_edge
	tot0 = tot_edge2
	if len(G2.nodes()) > len(G1.nodes()):
		tot = tot_edge2
		tot0 = tot_edge
	c3 = len(G12.edges()) * 1.0 / tot
	c4 = len(G12.edges()) * 1.0 / tot0
	if len(G1.nodes()) < len(G2.nodes()):
		c10 = len(G12.edges()) * 1.0 / len(G1.edges())
		EGf = len(G2.subgraph(nodes2).edges())
		c11 = len(G12.edges()) * 1.0 / (len(G1.edges()) + EGf - len(G12.edges()) * 1.0)
	else:
		c10 = len(G12.edges()) * 1.0 / len(G2.edges())
		EGf = len(G1.subgraph(nodes).edges())
		c11 = len(G12.edges()) * 1.0 / (len(G2.edges()) + EGf - len(G12.edges()) * 1.0)
	return [c10, c3, c11]

def intersec(G1, G2, matching):
	intersec = nx.Graph()
	for e in G1.edges():
		e0 = e[0]
		e1 = e[1]
		if e0 in matching and e1 in matching:
			if G2.has_edge(matching[e0], matching[e1]):
				intersec.add_edge(e0, e1)
	return intersec

def comms_has_edges(G, comm1, comm2):
	for node1 in comm1:
		for node2 in comm2: 
			if G.has_edge(node1, node2):
				return True
	return False


def EC(G1, G2, matching):
	intersec = nx.Graph()
	ec = 0
	for e in G1.edges():
		e0 = e[0]
		e1 = e[1]
		if e0 in matching and e1 in matching:
			if G2.has_edge(matching[e0], matching[e1]):
				ec += 1
	return ec / (len(G1.edges()) * 1.0), ec

def sortdict(d):
	s = d.items()
	s.sort()
	print '{' + ', '.join(map(lambda t: ': '.join(map(repr, t)), s)) + '}'


def load_data3(sp1, sp2, sp3):
	filename = '/home/90days/ekazemi/blast/proc/'
	dict1_prot_c = pickle.load(open(filename + sp1+'_prot_to_c.dat', 'rb'))
	dict1_c_prot = pickle.load(open(filename + sp1+'_c_to_prot.dat', 'rb'))
	dict2_prot_c = pickle.load(open(filename + sp2+'_prot_to_c.dat', 'rb'))
	dict2_c_prot = pickle.load(open(filename + sp2+'_c_to_prot.dat', 'rb'))
	dict3_prot_c = pickle.load(open(filename + sp3+'_prot_to_c.dat', 'rb'))
	dict3_c_prot = pickle.load(open(filename + sp3+'_c_to_prot.dat', 'rb'))
	
	G1 = load_graph(sp1)
	G2 = load_graph(sp2)
	G3 = load_graph(sp3)
	return G1, G2, G3, dict1_prot_c, dict1_c_prot, dict2_prot_c, dict2_c_prot, dict3_prot_c, dict3_c_prot

def load_graph(sp):
	f = open( '../proc/'+sp+'.txt', 'r' )
	G = nx.Graph()
	cnt = 0
	for line in f:
		if line.find( '#' ) != -1: continue
		cnt+=1
		if cnt == 1: continue
		parsed = line.replace( "\n", "" ).split()
		if parsed[0] != parsed[1]:
			G.add_edge(int(parsed[0]),int(parsed[1]))
	return G

def read_blast(sp1, sp2, min_blast):
	f = open( '/Users/ekazemi/Desktop/sorted_blast_' + sp1 +'-' + sp2 +'.blast', 'r' )
	lr = defaultdict(list)
	rl = defaultdict(list)
	pairs_dict12 = dict()
	pairs_dict21 = dict()
	for line in f:
		parsed = line.replace( "\n", "" ).split()
		if int(parsed[6]) >= min_blast:
			prot1 = parsed[1]
			prot2 = parsed[4]
			pairs_dict12[(prot1, prot2)] = int(parsed[6])
			pairs_dict21[(prot2, prot1)] = int(parsed[6])
			lr[prot1].append((prot2, int(parsed[6])))
			rl[prot2].append((prot1, int(parsed[6])))	
		else:
			break
	print len(lr), min_blast
	return lr, rl, pairs_dict12, pairs_dict21

def load_data4(sp1, sp2, sp3, sp4):
	filename = '/home/90days/ekazemi/blast/proc/'
	dict1_prot_c = pickle.load(open(filename + sp1+'_prot_to_c.dat', 'rb'))
	dict1_c_prot = pickle.load(open(filename + sp1+'_c_to_prot.dat', 'rb'))
	dict2_prot_c = pickle.load(open(filename + sp2+'_prot_to_c.dat', 'rb'))
	dict2_c_prot = pickle.load(open(filename + sp2+'_c_to_prot.dat', 'rb'))
	dict3_prot_c = pickle.load(open(filename + sp3+'_prot_to_c.dat', 'rb'))
	dict3_c_prot = pickle.load(open(filename + sp3+'_c_to_prot.dat', 'rb'))
	dict4_prot_c = pickle.load(open(filename + sp4+'_prot_to_c.dat', 'rb'))
	dict4_c_prot = pickle.load(open(filename + sp4+'_c_to_prot.dat', 'rb'))
	
	G1 = load_graph(sp1)
	G2 = load_graph(sp2)
	G3 = load_graph(sp3)
	G4 = load_graph(sp4)

	return G1, G2, G3, G4, dict1_prot_c, dict1_c_prot, dict2_prot_c, dict2_c_prot, dict3_prot_c, dict3_c_prot, dict4_prot_c, dict4_c_prot

def load_data5(sp1, sp2, sp3, sp4, sp5):
	filename = '/home/90days/ekazemi/blast/proc/'
	dict1_prot_c = pickle.load(open(filename + sp1+'_prot_to_c.dat', 'rb'))
	dict1_c_prot = pickle.load(open(filename + sp1+'_c_to_prot.dat', 'rb'))
	dict2_prot_c = pickle.load(open(filename + sp2+'_prot_to_c.dat', 'rb'))
	dict2_c_prot = pickle.load(open(filename + sp2+'_c_to_prot.dat', 'rb'))
	dict3_prot_c = pickle.load(open(filename + sp3+'_prot_to_c.dat', 'rb'))
	dict3_c_prot = pickle.load(open(filename + sp3+'_c_to_prot.dat', 'rb'))
	dict4_prot_c = pickle.load(open(filename + sp4+'_prot_to_c.dat', 'rb'))
	dict4_c_prot = pickle.load(open(filename + sp4+'_c_to_prot.dat', 'rb'))
	dict5_prot_c = pickle.load(open(filename + sp5+'_prot_to_c.dat', 'rb'))
	dict5_c_prot = pickle.load(open(filename + sp5+'_c_to_prot.dat', 'rb'))

	G1 = load_graph(sp1)
	G2 = load_graph(sp2)
	G3 = load_graph(sp3)
	G4 = load_graph(sp4)
	G5 = load_graph(sp5)

	return G1, G2, G3, G4, G5, dict1_prot_c, dict1_c_prot, dict2_prot_c, dict2_c_prot, dict3_prot_c, dict3_c_prot, dict4_prot_c, dict4_c_prot, dict5_prot_c, dict5_c_prot

def	find_conserved_3(G1, G2, G3, match12, match13):
	intersec = nx.Graph()
	for e in G1.edges():
		if e[0]in match12 and e[0] in match13 and e[1]in match12 and e[1] in match13:
			e02 = match12[e[0]]
			e03 = match13[e[0]]
			e12 = match12[e[1]]
			e13 = match13[e[1]]
			if G2.has_edge(e02, e12) and G2.has_edge(e02, e12):
				intersec.add_edge(e[0], e[1])
	return intersec


def load_graph_int(sp1):
	filename = '/home/90days/ekazemi/blast/proc/'
	dict1_prot_c = pickle.load(open(filename + sp1+'_prot_to_c.dat', 'rb'))
	dict1_c_prot = pickle.load(open(filename + sp1+'_c_to_prot.dat', 'rb'))
	
	f = open( filename + sp1+'.txt', 'r' )
	G1 = nx.Graph()
	cnt = 0
	for line in f:
		if cnt != 0:	
			parsed = line.replace( "\n", "" ).split()
			if parsed[0] != parsed[1]:
				G1.add_edge(int(parsed[0]),int(parsed[1]))
		cnt+=1
	
	return G1, dict1_prot_c, dict1_c_prot

def load_go(name):
	filename = '/home/90days/ekazemi/blast/proc/'
	all_go = defaultdict(list)
	sp = 'ce'
	GO = pickle.load(open(sp +'_prot_'+ name + 'go.pkl', 'rb'))
		
	dict_prot_c = pickle.load(open(filename + sp+'_prot_to_c.dat', 'rb'))
	for prot in GO:
		all_go[sp + str(dict_prot_c[prot])] = GO[prot]

	sp = 'dm'
	GO = pickle.load(open(sp +'_prot_'+ name + 'go.pkl', 'rb'))
	dict_prot_c = pickle.load(open(filename + sp+'_prot_to_c.dat', 'rb'))
	for prot in GO:
		all_go[sp + str(dict_prot_c[prot])] = GO[prot]

	sp = 'hs'
	GO = pickle.load(open(sp +'_prot_'+ name + 'go.pkl', 'rb'))
	dict_prot_c = pickle.load(open(filename + sp+'_prot_to_c.dat', 'rb'))
	for prot in GO:
		all_go[sp + str(dict_prot_c[prot])] = GO[prot]

	sp = 'mm'
	GO = pickle.load(open(sp +'_prot_'+ name + 'go.pkl', 'rb'))
	dict_prot_c = pickle.load(open(filename + sp+'_prot_to_c.dat', 'rb'))
	for prot in GO:
		all_go[sp + str(dict_prot_c[prot])] = GO[prot]

	sp = 'sc'
	GO = pickle.load(open(sp +'_prot_'+ name + 'go.pkl', 'rb'))
	dict_prot_c = pickle.load(open(filename + sp+'_prot_to_c.dat', 'rb'))
	for prot in GO:
		all_go[sp + str(dict_prot_c[prot])] = GO[prot]

	pickle.dump(all_go, open(name + "go.pkl", "wb"))

	return pickle.load(open(name + "go.pkl", "rb"))


def mean_entropy(cluster, all_go):
	prots = set()
	go_count = dict()
	for prot in cluster:
		if prot in all_go:
			prots.add(prot)
			go_prots = set(all_go[prot])
			for go in go_prots:
				if go not in go_count:
					go_count[go] = 1
				else:
					go_count[go] = go_count[go] + 1
					#print 'hi'
	if len(prots) > 1:
		size = len(prots) * 1.0
		entropy = 0
		d = len(go_count)
		for go in go_count:
			cnt = go_count[go]
			p = cnt / size
			#print p
			entropy += (-1.0 * p * math.log(p))

		if d > 1:
			nentropy = entropy / math.log(d)
			#nentropy = entropy / d
		else:
			nentropy = 0
		entropy = entropy / math.log(2)
		#nentropy = nentropy / math.log(2)
		return entropy, nentropy
	else:
		return -1, -1

def write_go_dict(filename, name):
	sps = ['ce', 'dm', 'hs', 'mm', 'sc']
	all_go = defaultdict(set)
	f = open(filename, 'r')
	for line in f:
		parsed = line.replace('\n', '').split('\t')
		for i in range(1, len(parsed)):
			if parsed[i].find('GO') != -1:
				all_go[parsed[0]].add(parsed[i])
	f.close()
	for sp in sps:
		graph_file = '../proc/' + sp + '.interaction'
		prots = set()
		f = open(graph_file, 'r')
		for line in f:
			parsed = line.replace('\n', '').split('\t')
			prots.add(parsed[0])
			prots.add(parsed[1])
		cnt = 0
		prot_go = defaultdict(list)
		for prot in prots:
			if prot in all_go:
				cnt += 1
				for go in all_go[prot]:
					if go.find('GO') != -1:
						prot_go[prot].append(go)
		gofile = sp + '_prot_' + name + 'go.pkl'
		pickle.dump(prot_go, open(gofile, 'wb'))
		print cnt	
		


def four(sp_seq, lr, nodes):
	quints = []
	prots = ['-1', '-1', '-1', '-1', '-1']
	prot1 = '-1'
	prot2 = '-1'
	prot3 = '-1'
	prot4 = '-1'
	prot5 = '-1'
	cnt = 0
	for prot1 in lr[sp_seq[0]][sp_seq[1]]:
		flag = False
		#print prot1, counter, size
		if prot1 not in nodes[sp_seq[0]]:
			for prot_val2 in lr[sp_seq[0]][sp_seq[1]][prot1]:
				if flag:
					break;
				prot2 = prot_val2[0]
				val12 = prot_val2[1]
				if prot2 in lr[sp_seq[1]][sp_seq[2]]:
					if prot2 not in nodes[sp_seq[1]]:
						for prot_val3 in lr[sp_seq[1]][sp_seq[2]][prot2]:
							if flag:
								break;
							prot3 = prot_val3[0]
							val23 = prot_val3[1]
							if prot3 in lr[sp_seq[2]][sp_seq[3]]:
								if prot3 not in nodes[sp_seq[2]]:
									for prot_val4 in lr[sp_seq[2]][sp_seq[3]][prot3]:
										prot4 = prot_val4[0]
										val34 = prot_val4[1]
										if prot4 not in nodes[sp_seq[3]]:
											cnt += 1
											val = val12 + val23 + val34
											nodes[sp_seq[0]].add(prot1)
											nodes[sp_seq[1]].add(prot2)
											nodes[sp_seq[2]].add(prot3)
											nodes[sp_seq[3]].add(prot4)
											prots[sp_seq[0]] = prot1
											prots[sp_seq[1]] = prot2
											prots[sp_seq[2]] = prot3
											prots[sp_seq[3]] = prot4
											quints.append([prots[0], prots[1], prots[2], prots[3], prots[4], val])
											flag = True
											break;
	print 'start sorting!', len(quints), 4
	return quints, nodes
def three(sp_seq, lr, nodes):
	quints = []
	prot1 = '-1'
	prot2 = '-1'
	prot3 = '-1'
	prot4 = '-1'
	prot5 = '-1'
	prots = ['-1', '-1', '-1', '-1', '-1']
	for prot1 in lr[sp_seq[0]][sp_seq[1]]:
		flag = False
		#print prot1, counter, size
		
		if prot1 not in nodes[sp_seq[0]]:
			for prot_val2 in lr[sp_seq[0]][sp_seq[1]][prot1]:
				if flag:
					break;
				prot2 = prot_val2[0]
				val12 = prot_val2[1]
				if prot2 in lr[sp_seq[1]][sp_seq[2]]:
					if prot2 not in nodes[sp_seq[1]]:
						for prot_val3 in lr[sp_seq[1]][sp_seq[2]][prot2]:
							if flag:
								break;
							prot3 = prot_val3[0]
							val23 = prot_val3[1]
							if prot3 not in nodes[sp_seq[2]]:
								val = val12 + val23
								nodes[sp_seq[0]].add(prot1)
								nodes[sp_seq[1]].add(prot2)
								nodes[sp_seq[2]].add(prot3)
								prots[sp_seq[0]] = prot1
								prots[sp_seq[1]] = prot2
								prots[sp_seq[2]] = prot3
								quints.append([prots[0], prots[1], prots[2], prots[3], prots[4], val])
								flag = True
								break;
	print 'start sorting!', len(quints), 3
	return quints, nodes

def two(sp_seq, lr, nodes):
	quints = []
	prot1 = '-1'
	prot2 = '-1'
	prot3 = '-1'
	prot4 = '-1'
	prot5 = '-1'
	prots = ['-1', '-1', '-1', '-1', '-1']
	for prot1 in lr[sp_seq[0]][sp_seq[1]]:
		flag = False
		#print prot1, counter, size
		if prot1 not in nodes[sp_seq[0]]:
			for prot_val2 in lr[sp_seq[0]][sp_seq[1]][prot1]:
				if flag:
					break;
				prot2 = prot_val2[0]
				val12 = prot_val2[1]
				if prot2 not in nodes[sp_seq[1]]:
					val = val12
					nodes[sp_seq[0]].add(prot1)
					nodes[sp_seq[1]].add(prot2)
					prots[sp_seq[0]] = prot1
					prots[sp_seq[1]] = prot2
					quints.append([prots[0], prots[1], prots[2], prots[3], prots[4], val])
					flag = True
					break;
	print 'start sorting!', len(quints), 2
	return quints, nodes


def five(sp_seq, lr, nodes):
	quints = []
	prot1 = '-1'
	prot2 = '-1'
	prot3 = '-1'
	prot4 = '-1'
	prot5 = '-1'
	prots = ['-1', '-1', '-1', '-1', '-1']
	cnt = 0
	tot = 0
	print len(lr[sp_seq[0]][sp_seq[1]])
	print len(lr[sp_seq[1]][sp_seq[2]])
	print len(lr[sp_seq[2]][sp_seq[3]])
	print len(lr[sp_seq[3]][sp_seq[4]])
	for prot1 in lr[sp_seq[0]][sp_seq[1]].keys():
		tot += 1
		flag = False
		if prot1 not in nodes[sp_seq[0]]:
			for prot_val2 in lr[sp_seq[0]][sp_seq[1]][prot1]:
				if flag:
					break;
				prot2 = prot_val2[0]
				val12 = prot_val2[1]
				if prot2 in lr[sp_seq[1]][sp_seq[2]]:
					if prot2 not in nodes[sp_seq[1]]:
						for prot_val3 in lr[sp_seq[1]][sp_seq[2]][prot2]:
							if flag:
								break;
							prot3 = prot_val3[0]
							val23 = prot_val3[1]
							if prot3 in lr[sp_seq[2]][sp_seq[3]]:
								if prot3 not in nodes[sp_seq[2]]:
									for prot_val4 in lr[sp_seq[2]][sp_seq[3]][prot3]:
										if flag:
											break;
										prot4 = prot_val4[0]
										val34 = prot_val4[1]
										if prot4 in lr[sp_seq[3]][sp_seq[4]]:
											if prot4 not in nodes[sp_seq[3]]:
												for prot_val5 in lr[sp_seq[3]][sp_seq[4]][prot4]:
													prot5 = prot_val5[0]
													val45 = prot_val5[1]
													if prot5 not in nodes[sp_seq[4]]:
														val = val12 + val23 + val34 + val45
														cnt += 1
														print val12, val23, val34, val45, cnt, tot
														nodes[sp_seq[0]].add(prot1)
														nodes[sp_seq[1]].add(prot2)
														nodes[sp_seq[2]].add(prot3)
														nodes[sp_seq[3]].add(prot4)
														nodes[sp_seq[4]].add(prot5)
														prots[sp_seq[0]] = prot1
														prots[sp_seq[1]] = prot2
														prots[sp_seq[2]] = prot3
														prots[sp_seq[3]] = prot4
														prots[sp_seq[4]] = prot5
														del lr[sp_seq[0]][sp_seq[1]][prot1]
														del lr[sp_seq[1]][sp_seq[2]][prot2]
														del lr[sp_seq[2]][sp_seq[3]][prot3]
														del lr[sp_seq[3]][sp_seq[4]][prot4]
														quints.append([prots[0], prots[1], prots[2], prots[3], prots[4], val])
														flag = True
														break;
													else:
														lr[sp_seq[3]][sp_seq[4]][prot4].remove((prot5, val45))
											else:
												del lr[sp_seq[3]][sp_seq[4]][prot4]
												lr[sp_seq[2]][sp_seq[3]][prot3].remove((prot4, val34))
								else:
									del lr[sp_seq[2]][sp_seq[3]][prot3]
									lr[sp_seq[1]][sp_seq[2]][prot2].remove((prot3, val23))
					else:
						del lr[sp_seq[1]][sp_seq[2]][prot2]
						lr[sp_seq[0]][sp_seq[1]][prot1].remove((prot2, val12))
		else:
			del lr[sp_seq[0]][sp_seq[1]][prot1]
	print len(lr[sp_seq[0]][sp_seq[1]])
	print len(lr[sp_seq[1]][sp_seq[2]])
	print len(lr[sp_seq[2]][sp_seq[3]])
	print len(lr[sp_seq[3]][sp_seq[4]])
	print 'start sorting!', len(quints), 5
	return quints, nodes, lr



def load_fuse_graph_int(sp1):
	filename = '/home/90days/ekazemi/blast/proc/'
	dict1_prot_c = pickle.load(open(filename + sp1+'_prot_to_c.dat', 'rb'))
	dict1_c_prot = pickle.load(open(filename + sp1+'_c_to_prot.dat', 'rb'))
	
	f = open( filename + sp1 + '.txt', 'r' )
	G1 = nx.Graph()
	cnt = 0
	for line in f:
		if cnt != 0:	
			parsed = line.replace( "\n", "" ).split()
			if parsed[0] != parsed[1]:
				G1.add_edge(int(parsed[0]),int(parsed[1]))
		cnt+=1
	
	return G1, dict1_prot_c, dict1_c_prot


def read_fuse_blast(sp1, sp2, min_blast):
	f = open( '../proc/'+ sp1 +'_' + sp2 +'.sim', 'r' )
	lr = defaultdict(list)
	rl = defaultdict(list)
	pairs_dict12 = dict()
	pairs_dict21 = dict()
	for line in f:
		parsed = line.replace( "\n", "" ).split('\t')
		if float(parsed[2]) >= min_blast:
			prot1 = parsed[0]
			prot2 = parsed[1]
			pairs_dict12[(prot1, prot2)] = float(parsed[2])
			pairs_dict21[(prot2, prot1)] = float(parsed[2])
			lr[prot1].append((prot2, float(parsed[2])))
			rl[prot2].append((prot1, float(parsed[2])))	
		
	print len(lr), min_blast
	return lr, rl, pairs_dict12, pairs_dict21 


def CIQ(clusters, G1, G2, G3, G4, G5):
	len_clusters = len(clusters)
	val = 0
	for i in range(len_clusters):
		for j in range(i+1, len_clusters):
			cluster1 = clusters[i]
			cluster2 = clusters[j]
			snm = 0
			for index in range(5):
				if cluster1[index] != -1 and cluster2[index] != -1:
					snm+=1
			
			if snm > 1:
				snm_ = 0
				index = 0
				if cluster1[index] != -1 and cluster2[index] != -1:
					if G1.has_edge(cluster1[index], cluster2[index]):
						snm_ += 1.
				index = 1
				if cluster1[index] != -1 and cluster2[index] != -1:
					if G2.has_edge(cluster1[index], cluster2[index]):
						snm_ += 1.
				index = 2
				if cluster1[index] != -1 and cluster2[index] != -1:
					if G3.has_edge(cluster1[index], cluster2[index]):
						snm_ += 1.
				index = 3
				if cluster1[index] != -1 and cluster2[index] != -1:
					if G4.has_edge(cluster1[index], cluster2[index]):
						snm_ += 1.
				index = 4
				if cluster1[index] != -1 and cluster2[index] != -1:
					if G5.has_edge(cluster1[index], cluster2[index]):
						snm_ += 1.
				if snm_ > 1:
					val += (snm_ / snm)
		#if snm != 0 and i % 100 == 0:
			#print snm, snm_, val, i, val / ((i+1) * len_clusters)

def CIQ_general(clusters, G1, G2, G3, G4, G5):
	garph_dict = dict()
	garph_dict['ce'] = G1
	garph_dict['dm'] = G2
	garph_dict['hs'] = G3
	garph_dict['mm'] = G4
	garph_dict['sc'] = G5

	len_clusters = len(clusters)

	val = 0
	CIQ_n = 0
	modified_clusters = dict()
	set_clusters = defaultdict(set)
	for i in range(len_clusters):
		modified_clusters[i] = defaultdict(list)
		for prot in clusters[i]:
				modified_clusters[i][prot[0:2]].append(prot)
				set_clusters[i].add(prot[0:2])

	enm_tot = 0
	for i in range(len_clusters):
		for j in range(i+1, len_clusters):
			intersec = ( set_clusters[i] & set_clusters[j])
			snm = len(intersec)
			#if snm > 1:
			snm_ = 0
			enm = 0
			for sp in intersec:
				flag = False
				for prot1 in modified_clusters[i][sp]:
					for prot2 in modified_clusters[j][sp]:
						prot1_int = int(prot1[2:])
						prot2_int = int(prot2[2:])
						if garph_dict[sp].has_edge(prot1_int, prot2_int):
							if not flag:
								snm_ += 1.0
								flag = True
							enm += 1.0
			if snm > 1:
				if snm_ > 1:
					cnm = snm_ / snm
					val += cnm
					CIQ_n += (cnm * enm)
			enm_tot += enm
					#print modified_clusters[i][sp], modified_clusters[j][sp]
		#if snm != 0 and i % 100 == 0:
			#print snm, snm_, val, i, val / ((i+1) * len_clusters), CIQ_n / enm_tot
	#print '++++++++++++++++++', CIQ_n / enm_tot
	return CIQ_n / enm_tot
	


def blast(sp1, sp2, node1, node2, out):
	species = {'sc':'cerevisiae','ce':'elegans','hs':'homosapiens','dm':'melanogaster','mm':'musculus'}
	FNULL = open(os.devnull, 'w')
	retcode = subprocess.call(["./blastp", "-query", "../data/"+species[sp1]+"_fasta/"+node1+".txt", "-subject", "../data/"+species[sp2]+"_fasta/"+node2+".txt", "-out" , ""+out+".txt"],stdout=FNULL, stderr=subprocess.STDOUT)
	ft= open( ""+out+".txt", "r" )
	s_flag = True
	val = -1
	for line in ft:
		if not s_flag:
			break
		if line.find("Length=") != -1:
			line = line.replace( "\n", "" )
		if line.find("Score = ") != -1:
			line = line.replace( "\n", "" )
			line = line[9:]
			ind_bits = line.index("b")
			val =  float(line[:ind_bits-1].replace( " ", "" ))
			s_flag = False	
	return 	val
	ft.close()


def load_all_data(min_blast, sps):
	vector = []
	for i in range(len(sps)):
		vector.append([])
	lr = []
	rl = []
	pairs_dict = []
	nodes = []
	for i in range(len(sps)):
		nodes.append(set())
		lr.append(copy.copy(vector))
		rl.append(copy.copy(vector))
		pairs_dict.append(copy.copy(vector))

	for i in range(len(sps)):
		for j in range(i+1, len(sps)):
			lrij, rlij, pairs_dictij12, pairs_dictij21= read_blast(sps[i], sps[j], min_blast)
			lr[i][j] = lrij
			rl[i][j] = rlij
			lr[j][i] = rlij
			rl[j][i] = lrij
			pairs_dict[i][j] = pairs_dictij12
			pairs_dict[j][i] = pairs_dictij21
	return lr, rl, pairs_dict


def load_all_prot_seq(sps):
	species = {'sc':'cerevisiae','ce':'elegans','hs':'homosapiens','dm':'melanogaster','mm':'musculus'}
	all_prot_seq = []
	for sp in sps:
		all_prot = set()
		f = open( "../data/prot_names/"+species[sp]+".pr", "r" )
		for line in f:
			prot = line.replace( "\n", "" )
			all_prot.add(prot)
		f.close()
		all_prot_seq.append(all_prot)
	return all_prot_seq


def bipartiteMatch(graph):
	'''Find maximum cardinality matching of a bipartite graph (U,V,E).
	The input format is a dictionary mapping members of U to a list
	of their neighbors in V.  The output is a triple (M,A,B) where M is a
	dictionary mapping members of V to their matches in U, A is the part
	of the maximum independent set in U, and B is the part of the MIS in V.
	The same object may occur in both U and V, and is treated as two
	distinct vertices if this happens.'''
	
	# initialize greedy matching (redundant, but faster than full search)
	matching = {}
	for u in graph:
		for v in graph[u]:
			if v not in matching:
				matching[v] = u
				break
	
	while 1:
		# structure residual graph into layers
		# pred[u] gives the neighbor in the previous layer for u in U
		# preds[v] gives a list of neighbors in the previous layer for v in V
		# unmatched gives a list of unmatched vertices in final layer of V,
		# and is also used as a flag value for pred[u] when u is in the first layer
		preds = {}
		unmatched = []
		pred = dict([(u,unmatched) for u in graph])
		for v in matching:
			del pred[matching[v]]
		layer = list(pred)
		
		# repeatedly extend layering structure by another pair of layers
		while layer and not unmatched:
			newLayer = {}
			for u in layer:
				for v in graph[u]:
					if v not in preds:
						newLayer.setdefault(v,[]).append(u)
			layer = []
			for v in newLayer:
				preds[v] = newLayer[v]
				if v in matching:
					layer.append(matching[v])
					pred[matching[v]] = v
				else:
					unmatched.append(v)
		
		# did we finish layering without finding any alternating paths?
		if not unmatched:
			unlayered = {}
			for u in graph:
				for v in graph[u]:
					if v not in preds:
						unlayered[v] = None
			return (matching,list(pred),list(unlayered))

		# recursively search backward through layers to find alternating paths
		# recursion returns true if found path, false otherwise
		def recurse(v):
			if v in preds:
				L = preds[v]
				del preds[v]
				for u in L:
					if u in pred:
						pu = pred[u]
						del pred[u]
						if pu is unmatched or recurse(pu):
							matching[v] = u
							return 1
			return 0

		for v in unmatched: recurse(v)

	return matching


def read_DAG():
	filename = 'go-basic.obo'
	#filename = 'go.obo'
	DG=nx.DiGraph()

	f = open(filename, 'r')
	term_flag = False
	id_flag = False
	cnt = 0 
	for line in f:
		cnt += 1
		proc = line.replace('\n', '')
		if proc.find('[Term]') != -1:
			term_flag = True
			id_flag = False

		if term_flag:
			if proc.find('id: GO:') != -1:
				term_flag = False
				id_flag = True
				go_id = proc[4:14]
				DG.add_node(go_id)
		else:
			if id_flag:
				if proc.find('is_a: GO:') != -1:
					go_parent = proc[6:16]
					#print go_id, go_parent
					DG.add_edge(go_parent,go_id)
	print nx.is_directed_acyclic_graph(DG)
	return DG

def go_level(DG):
	roots = set()
	go_level = dict()
	level = 1
	for node in DG.nodes():
		if len(DG.predecessors(node)) == 0 and len(DG.successors(node)) != 0:
			roots.add(node)
			go_level[node] = level
		if len(DG.predecessors(node)) == 0 and len(DG.successors(node)) == 0:
			go_level[node] = level
	while(len(roots)) > 0:
		level += 1
		roots_tmp = set()
		for node in roots:
			for go_child in DG.successors(node):
				if go_child not in go_level:
					go_level[go_child] = level
					roots_tmp.add(go_child)
		roots = copy.copy(roots_tmp)
	return go_level

def go_IC():
	"""name = 'AD'
	PgoA = pickle.load(open("GO"+ name +".dat", "rb"))
	goA = set()
	for prot in PgoA:
		for go in PgoA[prot]:
			goA.add(go)
	nb_A = len(PgoA) * 1.0

	name = 'PD'
	PgoP = pickle.load(open("GO"+ name +".dat", "rb"))
	goP = set()
	for prot in PgoP:
		for go in PgoP[prot]:
			goP.add(go)
	nb_P = len(PgoP) * 1.0

	name = 'FD'
	PgoF = pickle.load(open("GO"+ name +".dat", "rb"))
	goF = set()
	for prot in PgoF:
		for go in PgoF[prot]:
			goF.add(go)
	nb_F = len(PgoF) * 1.0

	name = 'CD'
	PgoC = pickle.load(open("GO"+ name +".dat", "rb"))
	goC = set()
	for prot in PgoC:
		for go in PgoC[prot]:
			goC.add(go)
	nb_C = len(PgoC) * 1.0
	

	go_cnt = dict()
	for prot in PgoA:
		for go in PgoA[prot]:
			if go != 'GO:0008150' and go != 'GO:0005575' and go != 'GO:0003674': 
				if go not in go_cnt:
					go_cnt[go] = 1
				else:
					go_cnt[go] += 1
	go_IC = dict()
	cnt = 0
	val = 0
	for go in go_cnt:
		val += 1
		if go in goP:
			go_IC[go] = -1.0 * math.log(go_cnt[go] / nb_P) / math.log(2)
		else:
			if go in goF:
				go_IC[go] = -1.0 * math.log(go_cnt[go] / nb_F) / math.log(2)
			else:
				if go in goC:
					go_IC[go] = -1.0 * math.log(go_cnt[go] / nb_C) / math.log(2)
				else:
					cnt +=1
					print 'Shit!', cnt, val
					exit()
	print nb_A, nb_C, nb_F, nb_P
	pickle.dump(go_IC, open('go_IC.dat', 'wb'))"""
	go_IC = pickle.load(open('IC.dat', 'rb'))

	return go_IC
		
def read_alignments(alg,min_blast):
	sps = ['ce', 'dm', 'hs', 'mm', 'sc']

	if alg.find('.txt') != -1:
		matchfile = alg
		f = open(matchfile, 'r')
		all_clusters = []
		clusters = []
		for line in f:
			parsed = line.replace('\n', '').split('\t')
			values = []
			for val in parsed:
				values.append(int(val))
			#prot1 = int(parsed[0])
			#prot2 = int(parsed[1])
			#prot3 = int(parsed[2])
			#prot4 = int(parsed[3])
			#prot5 = int(parsed[4])
			#values = [prot1, prot2, prot3, prot4, prot5]
			cluster = []
			for i in range(len(values)):
				if values[i] != -1:
					cluster.append(sps[i] + str(values[i]))
			#print cluster
			all_clusters.append(cluster)
		return all_clusters

	if alg == '1' or alg == '0' or alg == 'pgm' or alg == 'pgmc':
		matchfile = 'match-'+alg+'-'+sps[0]+'-'+sps[1]+'-'+sps[2]+'-'+sps[3]+'-'+sps[4]+'-'+str(min_blast)+'.txt'
		f = open(matchfile, 'r')
		cnt = 0
		all_clusters = []
		clusters = []
		for line in f:
			parsed = line.replace('\n', '').split('\t')
			prot1 = int(parsed[0])
			prot2 = int(parsed[1])
			prot3 = int(parsed[2])
			prot4 = int(parsed[3])
			prot5 = int(parsed[4])
			cnt+=1
			values = [prot1, prot2, prot3, prot4, prot5]	
			cluster = []
			for i in range(len(sps)):
				if values[i] != -1:
					cluster.append(sps[i] + str(values[i]))
			all_clusters.append(cluster)
		return all_clusters
	if alg == 'R':
		alg = 'pgm'
		matchfile = 'match-'+alg+'-'+sps[0]+'-'+sps[1]+'-'+sps[2]+'-'+sps[3]+'-'+sps[4]+'-'+str(min_blast)+'.txt'
		f = open(matchfile, 'r')
		cnt = 0
		all_clusters = []
		clusters = []
		for line in f:
			parsed = line.replace('\n', '').split('\t')
			prot1 = int(parsed[0])
			prot2 = int(parsed[1])
			prot3 = int(parsed[2])
			prot4 = int(parsed[3])
			prot5 = int(parsed[4])
			cnt+=1
			values = [prot1, prot2, prot3, prot4, prot5]	
			cluster = []
			cnt = 0
			for i in range(len(sps)):
				if values[i] != -1:
					cnt+=1
					cluster.append(sps[i] + str(values[i]))
			if cnt != 5 or random.uniform(0, 1) < 0.3:
				all_clusters.append(cluster)
			else:
				cluster1 = random.sample(cluster, 2)
				cluster2 = list(set(cluster) - set(cluster1))
				all_clusters.append(cluster1)
				all_clusters.append(cluster2)
		return all_clusters
	if alg == 'B':
		matchfile = '../BEAMS/many_0' + min_blast + '.txt'
		f = open(matchfile, 'r')
		cnt = 0
		counts = np.zeros(5)
		all_clusters = []
		nodes = set()
		clusters = []
		all_prots = set()
		new_prot = 0
		for line in f:
			parsed = line.replace('\n', '').split(' ')
			cluster = []
			for i in range(len(parsed) - 1):
				cluster.append(parsed[i])
			#print cluster, parsed
			all_clusters.append(cluster)
		return all_clusters

	if alg == 'BB':
		matchfile = '../BEAMS/out_0' + min_blast + '.txt'
		f = open(matchfile, 'r')
		cnt = 0
		counts = np.zeros(5)
		all_clusters = []
		nodes = set()
		clusters = []
		all_prots = set()
		new_prot = 0
		for line in f:
			parsed = line.replace('\n', '').split(' ')
			cluster = []
			for i in range(len(parsed) - 1):
				cluster.append(parsed[i])
			#print cluster, parsed
			all_clusters.append(cluster)
		return all_clusters

	if alg == 'F':
		matchfile = '../FUSE/DATABASE/an_1000_output_file'
		val0 = 0
		val1 = 4950
		val2 = 4950 + 8532
		val3 = 4950 + 8532 + 19141
		val4 = 4950 + 8532 + 19141 + 10765
		f = open(matchfile, 'r')
		cnt = 0
		counts = np.zeros(5)
		all_clusters = []
		clusters = []
		nodes = set()
		tot_prots = 0
		for line in f:
			prot1 = -1
			prot2 = -1
			prot3 = -1
			prot4 = -1
			prot5 = -1
			parsed = line.replace('\n', '').split('\t')
			for i in range(2, len(parsed) - 1):
				#print parsed[i],
				nodes.add(parsed[i])
				prot = int(parsed[i])
				if prot >= val0 and prot < val1:
					prot1 = prot - val0
				if prot >= val1 and prot < val2:
					prot2 = prot - val1
				if prot >= val2 and prot < val3:
					prot3 = prot - val2
				if prot >= val3 and prot < val4:
					prot4 = prot - val3
				if prot >= val4:
					prot5 = prot - val4
			#print prot1, prot2, prot3, prot4, prot5
			values = [prot1, prot2, prot3, prot4, prot5]
		
			cnt = 0
			for val in values:
				if val == -1:
					cnt+=1
				else:
					tot_prots += 1
			cluster = []
			prot = prot1
			sp = sps[0]
			if prot != -1:
				cluster.append(sp + str(prot))
			prot = prot2
			sp = sps[1]
			if prot != -1:
				cluster.append(sp + str(prot))
			prot = prot3
			sp = sps[2]
			if prot != -1:
				cluster.append(sp + str(prot))
			prot = prot4
			sp = sps[3]
			if prot != -1:
				cluster.append(sp + str(prot))
			prot = prot5
			sp = sps[4]
			if prot != -1:
				cluster.append(sp + str(prot))
			#print cluster
			all_clusters.append(cluster)
		return all_clusters

	if alg == 'S':
		matchfile = '../SMETANA/SMETANA_ce_dm_hs_mm_sc.txt'
		#print matchfile
		f = open(matchfile, 'r')
		all_prots = defaultdict(set)
		prots_cnt = dict()
		names = {'ce': 0, 'dm': 1, 'hs': 2, 'mm': 3, 'sc': 4}
		sp_cnt = dict()
		sp_cnt[0] = 0 
		sp_cnt[1] = 0
		sp_cnt[2] = 0
		sp_cnt[3] = 0
		sp_cnt[4] = 0
		sp_cnt[5] = 0
		all_clusters = []
		for line in f:
			values = np.zeros(5)
			parsed = line.replace('\n', '').split(' ')
			cluster = []
			if len(parsed) < 4 or True:
				for prot in parsed:
					sp = prot[0:2]
					if sp in names:
						prot_int = int(prot[2:]) - 1
						all_prots[sp].add(prot_int)
						values[names[sp]] += 1
						cluster.append(sp + str(prot_int))
				all_clusters.append(cluster)
		return all_clusters

	if alg == 'C':
		matchfile = '../SMETANA/CSRW_ce_dm_hs_mm_sc.txt'
		#print matchfile
		f = open(matchfile, 'r')
		all_prots = defaultdict(set)
		prots_cnt = dict()
		names = {'ce': 0, 'dm': 1, 'hs': 2, 'mm': 3, 'sc': 4}
		sp_cnt = dict()
		sp_cnt[0] = 0 
		sp_cnt[1] = 0
		sp_cnt[2] = 0
		sp_cnt[3] = 0
		sp_cnt[4] = 0
		sp_cnt[5] = 0
		all_clusters = []
		for line in f:
			values = np.zeros(5)
			parsed = line.replace('\n', '').split(' ')
			cluster = []
			if len(parsed) < 4 or True:
				for prot in parsed:
					sp = prot[0:2]
					if sp in names:
						prot_int = int(prot[2:]) - 1
						all_prots[sp].add(prot_int)
						values[names[sp]] += 1
						cluster.append(sp + str(prot_int))
				all_clusters.append(cluster)
		return all_clusters
	if alg == 'I':
		matchfile = '../iso/final_ce_dm_hs_mm_sc.txt'
		#print matchfile
		f = open(matchfile, 'r')
		all_prots = defaultdict(set)
		prots_cnt = dict()
		names = {'ce': 0, 'dm': 1, 'hs': 2, 'mm': 3, 'sc': 4}
		sp_cnt = dict()
		sp_cnt[0] = 0 
		sp_cnt[1] = 0
		sp_cnt[2] = 0
		sp_cnt[3] = 0
		sp_cnt[4] = 0
		sp_cnt[5] = 0
		all_clusters = []
		for line in f:
			values = np.zeros(5)
			parsed = line.replace('\n', '').split('\t')
			cluster = []
			if len(parsed) < 4 or True:
				for prot in parsed:
					sp = prot[0:2]
					if sp in names:
						prot_int = int(prot[2:]) - 1
						all_prots[sp].add(prot_int)
						values[names[sp]] += 1
						cluster.append(sp + str(prot_int))
				all_clusters.append(cluster)
		return all_clusters

def mna_ic(all_clusters,go_IC):
	sp1 = 'ce'
	sp2 = 'dm'
	sp3 = 'hs'
	sp4 = 'mm'
	sp5= 'sc'
	name = 'EAD'
	G1, dict1_prot_c, dict1_c_prot,= load_graph_int(sp1)
	net1 = set(G1.nodes())
	G2, dict2_prot_c, dict2_c_prot,= load_graph_int(sp2)
	net2 = set(G2.nodes())
	G3, dict3_prot_c, dict3_c_prot,= load_graph_int(sp3)
	net3 = set(G3.nodes())
	G4, dict4_prot_c, dict4_c_prot,= load_graph_int(sp4)
	net4 = set(G4.nodes())
	G5, dict5_prot_c, dict5_c_prot,= load_graph_int(sp5)

	all_go_prot_tmp = pickle.load(open("GO"+ name +".dat", "rb"))
	gos = defaultdict(set)
	prots = set()
	for prot in all_go_prot_tmp:
		for go in all_go_prot_tmp[prot]:
			gos[go].add(prot)
		prots.add(prot)



	nb_prot = len(prots)
	frequent_go = set()
	for go in gos:
		if  len(gos[go]) > nb_prot * 1.0:
			frequent_go.add(go)

	all_go_prot = defaultdict(set)
	for prot in all_go_prot_tmp:
		for go in all_go_prot_tmp[prot]:
			if go not in frequent_go:
				all_go_prot[prot].add(go)

	all_go = defaultdict(set)
	for prot in all_go_prot:
		prot_convert = '-1'
		if prot in dict1_prot_c:
			prot_convert = sp1 + str(dict1_prot_c[prot])
			for go in all_go_prot[prot]:
				all_go[prot_convert].add(go)
		if prot in dict2_prot_c:
			prot_convert = sp2 + str(dict2_prot_c[prot])
			for go in all_go_prot[prot]:
				all_go[prot_convert].add(go)
		if prot in dict3_prot_c:
			prot_convert = sp3 + str(dict3_prot_c[prot])
			for go in all_go_prot[prot]:
				all_go[prot_convert].add(go)
		if prot in dict4_prot_c:
			prot_convert = sp4 + str(dict4_prot_c[prot])
			for go in all_go_prot[prot]:
				all_go[prot_convert].add(go)
		if prot in dict5_prot_c:
			prot_convert = sp5 + str(dict5_prot_c[prot])
			for go in all_go_prot[prot]:
				all_go[prot_convert].add(go)
		if prot_convert == '-1':
			break


	all_prots = set()
	for cluster in all_clusters:
		for prot in cluster:
			all_prots.add(prot)

	tot_prots = len(all_prots)
	
	sp_int = {'ce': 0, 'dm': 1, 'hs': 2, 'mm': 3, 'sc': 4}
	tot_IC = 0
	IC = np.zeros(5)
	tot = 0
	ncnt_ent = 0
	counts = np.zeros(5)
	counts_annotated = np.zeros(5)
	counts_go = np.zeros(5)
	counts_annotated_go = np.zeros(5)
	ciq = 0
	correct_prots = 0
	annotated_prots = 0
	for cluster in all_clusters:
		
		inter_sec = []
		annotated = 0
		cluster_sp = np.zeros(5)
		cluster_sp_go = np.zeros(5)
		cluster_sp_cnt = np.zeros(5)
		cluster_sp_cnt_go = np.zeros(5)
		cnt = 0
		for prot in cluster:
			go = []
			cluster_sp[sp_int[prot[0:2]]] = 1
			cluster_sp_cnt[sp_int[prot[0:2]]] += 1
			if prot in all_go:
				cluster_sp_go[sp_int[prot[0:2]]] = 1
				cluster_sp_cnt_go[sp_int[prot[0:2]]] += 1
				go = all_go[prot]
				annotated += 1
				if cnt == 0:
					inter_sec = go
				else:
					inter_sec = intersect(go, inter_sec)
				cnt += 1
		
		if annotated > 1:				
			annotated_prots += len(cluster)
			if len(inter_sec) > 0:
				max_ic = 0
				for go in inter_sec:
					if go in go_IC:
						if go_IC[go] > max_ic:
							max_ic = go_IC[go]
				tot_IC += (max_ic * np.sum(cluster_sp))
				IC[5 - np.sum(cluster_sp)] += max_ic

			counts_annotated[5 - np.sum(cluster_sp)] += 1	
			tot+=1

	ic_str = (str(IC[3]) + '\t' + str(IC[2]) + '\t' + str(IC[1]) + '\t'+ str(IC[0]) + '\t' + str(counts_annotated[3]) + '\t' + str(counts_annotated[2]) + '\t' + str(counts_annotated[1]) + '\t'+ str(counts_annotated[0]))
	return ic_str	


def intersec_multiple(alg, min_blast):
	clusters = read_alignments(alg,min_blast)

	sps = ['ce', 'dm', 'hs', 'mm', 'sc']
	sp_int = {'ce': 0, 'dm': 1, 'hs': 2, 'mm': 3, 'sc': 4}

	G = [[], [], [], [], []]
	dict_c_prot = [[], [], [], [], []]
	G[0], dict1_prot_c, dict_c_prot[0]= load_graph_int(sps[0])
	G[1], dict2_prot_c, dict_c_prot[1]= load_graph_int(sps[1])
	G[2], dict3_prot_c, dict_c_prot[2]= load_graph_int(sps[2])
	G[3], dict4_prot_c, dict_c_prot[3]= load_graph_int(sps[3])
	G[4], dict5_prot_c, dict_c_prot[4]= load_graph_int(sps[4])

	cluster_separate = []
	for cluster in clusters:
		cluster_tmp = [[], [], [], [], []]
		for prot in cluster:
			sp = prot[0:2]
			ind = sp_int[sp]
			prot_int = int(prot[2:])
			cluster_tmp[ind].append(prot_int)
		cluster_separate.append(copy.copy(cluster_tmp))

	cluster_all_sp = []

	for cluster in cluster_separate:
		cnt = 0
		for i in range(5):
			if len(cluster[i]) > 0:
				cnt += 1
		if cnt == 5:
			cluster_all_sp.append(cluster)


	name = 'EAL5'
	all_go_prot = pickle.load(open("GO"+ name +".dat", "rb"))

	consistent_cluster = []
	cnt = 0
	tot_anot = 0
	c_anot = 0
	for cluster in clusters:
		for prot in cluster:
			sp = prot[0:2]
			ind = sp_int[sp]
			prot_int = int(prot[2:])
			if prot_int in dict_c_prot[ind]:
				if dict_c_prot[ind][prot_int] in all_go_prot:
					#print prot
					cnt +=1
				
		inter_sec = []
		annotated = 0
		cluster_sp = np.zeros(5)
		cnt = 0
		for prot in cluster:
			sp = prot[0:2]
			ind = sp_int[sp]
			prot_int = int(prot[2:])
			if prot_int in dict_c_prot[ind]:
				prot_name = dict_c_prot[ind][prot_int]
				go = []
				cluster_sp[ind] = 1
				if prot_name in all_go_prot:
					go = all_go_prot[prot_name]
					annotated += 1
					if cnt == 0:
						inter_sec = go
					else:
						inter_sec = intersect(go, inter_sec)
					cnt += 1
		
		if annotated > 1:
			tot_anot += 1
			if len(inter_sec) > 0:
				c_anot += 1
				cluster_tmp = [[], [], [], [], []]
				for prot in cluster:
					sp = prot[0:2]
					ind = sp_int[sp]
					prot_int = int(prot[2:])
					cluster_tmp[ind].append(prot_int)
				consistent_cluster.append(copy.copy(cluster_tmp))
	tot_edge = 0
	t_edge = 0
	conserved_edge = 0
	m_conserved_edge = 0
	cnt0 = 0
	edge_count = np.zeros(6)
	edge_count_tmp = np.zeros(6)
	aligned_edges = 0
	for ind1 in range(len(consistent_cluster)):
		for ind2 in range(ind1+1, len(consistent_cluster)):
			c1 = consistent_cluster[ind1]
			c2 = consistent_cluster[ind2]
			nb_edge = 0
			sp_cnt = 0
			sp_edge = 0
			for i in range(5):
				if len(c1[i]) > 0 and len(c2[i]) > 0:
					sp_cnt += 1
					edge_flag = False
					for p1 in c1[i]:
						for p2 in c2[i]:
							if G[i].has_edge(p1,p2):
								edge_flag = True
								nb_edge += 1
								if i == 0:
									#print i, p1, p2
									cnt0 += 1
					if edge_flag:
						sp_edge += 1
			if sp_cnt > 1:
				tot_edge += nb_edge
				if sp_cnt == sp_edge:
					t_edge += nb_edge
					conserved_edge += nb_edge
					edge_count[sp_cnt] += 1
				if nb_edge > 1:
					edge_count_tmp[sp_edge] += 1
					aligned_edges += nb_edge
				if sp_edge > 1:
					m_conserved_edge += nb_edge

	c_intersec = str(t_edge) + '\t' + str(tot_edge) + '\t' + str(m_conserved_edge) + '\t' + str(aligned_edges) + '\t' + str(conserved_edge)



	tot_edge = 0
	t_edge = 0
	conserved_edge = 0
	m_conserved_edge  = 0
	cnt0 = 0
	edge_count = np.zeros(6)
	edge_count_tmp = np.zeros(6)
	aligned_edges = 0
	for ind1 in range(len(cluster_separate)):
		for ind2 in range(ind1+1, len(cluster_separate)):
			c1 = cluster_separate[ind1]
			c2 = cluster_separate[ind2]
			nb_edge = 0
			sp_cnt = 0
			sp_edge = 0
			for i in range(5):
				if len(c1[i]) > 0 and len(c2[i]) > 0:
					sp_cnt += 1
					edge_flag = False
					for p1 in c1[i]:
						for p2 in c2[i]:
							if G[i].has_edge(p1,p2):
								edge_flag = True
								nb_edge += 1
								if i == 0:
									#print i, p1, p2
									cnt0 += 1
					if edge_flag:
						sp_edge += 1
			if sp_cnt > 1:
				tot_edge += nb_edge
				if sp_cnt == sp_edge:
					t_edge += nb_edge
					conserved_edge += nb_edge
					edge_count[sp_cnt] += 1
				if nb_edge > 1:
					edge_count_tmp[sp_edge] += 1
					aligned_edges += nb_edge
				if sp_edge > 1:
					m_conserved_edge += nb_edge
	t_intersec = str(t_edge) + '\t' + str(tot_edge) + '\t' + str(m_conserved_edge) + '\t' + str(aligned_edges) + '\t' + str(conserved_edge)
	return (c_intersec + '\t' + t_intersec)



def merge_clusters(clusters, val):

	shuffle(clusters)
	sps = ['ce', 'dm', 'hs', 'mm', 'sc']
	sp_int = {'ce': 0, 'dm': 1, 'hs': 2, 'mm': 3, 'sc': 4}

	G = [[], [], [], [], []]
	dict_c_prot = [[], [], [], [], []]
	G[0], dict1_prot_c, dict_c_prot[0]= load_graph_int(sps[0])
	G[1], dict2_prot_c, dict_c_prot[1]= load_graph_int(sps[1])
	G[2], dict3_prot_c, dict_c_prot[2]= load_graph_int(sps[2])
	G[3], dict4_prot_c, dict_c_prot[3]= load_graph_int(sps[3])
	G[4], dict5_prot_c, dict_c_prot[4]= load_graph_int(sps[4])


	matched = set()
	for c in clusters:
		for prot in c:
			matched.add(prot)
	unmatched = [[], [], [], [], []]
	for i in range(5):
		for node in G[i].nodes():
			prot = sps[i] + str(node)
			if prot not in matched:
				unmatched[i].append(prot)
	for i in range(5):
		shuffle(unmatched[i])
		print len(G[i].nodes()), len(unmatched[i])
	cnt = 0
	new_clusters = []
	for c in clusters:
		sp_inds = np.zeros(5)
		tmpc = copy.copy(c)
		for prot in c:
			sp = prot[0:2]
			sp_inds[sp_int[sp]] = 1
		for i in range(5):
			if sp_inds[i] == 0:
				if len(unmatched[i]) != 0:
					nprot = unmatched[i][0]
					unmatched[i].pop(0)
					tmpc.append(nprot)
		#print tmpc, c
		cnt +=1
		if cnt < val:
			new_clusters.append(tmpc)
		else:
			new_clusters.append(c)
	"""sps = ['ce', 'dm', 'hs', 'mm', 'sc']
	sp_int = {'ce': 0, 'dm': 1, 'hs': 2, 'mm': 3, 'sc': 4}
	lr, rl, pairs_dict = load_all_data(min_blast, sps)
	pairs = []
	for i in range(5):
		for j in range(i+1, 5):
			for pair in pairs_dict[i][j]:
				pairs.append([i, j, pair, pairs_dict[i][j][pair]])
	print len(pairs)
	pairs.sort(key = lambda tup: - 1 * tup[3])
	prot_ind = dict()
	nb_c = len(clusters)
	for i in range(nb_c):
		c1 = clusters[i]
		for prot in c1:
			prot_ind[prot] = i

	for pair in pairs:
		prot1 = sps[pair[0]] + pair[2][0]
		prot2 = sps[pair[1]] + pair[2][1]
		if prot1 in prot_ind and prot2 in prot_ind:
			if prot_ind[prot1] != prot_ind[prot2]:
				ind1 = prot_ind[prot1]
				ind2 = prot_ind[prot2]
				if ind1 < ind2:
					for prot in clusters[ind2]:
						prot_ind[prot] = ind1
				else:
					for prot in clusters[ind1]:
						prot_ind[prot] = ind2
	ind_int = dict()
	new_clusters = []
	cnt = 0
	
	for prot in prot_ind:
		ind = prot_ind[prot]
		if ind in ind_int:
			new_clusters[ind_int[ind]].append(prot)
		else:
			ind_int[ind] = cnt
			new_clusters.append([])
			new_clusters[cnt].append(prot)
			cnt+=1
	print len(new_clusters), len(clusters)
	clusters = []
	prot_ind = dict()
	cnt = 0
	for pair in pairs:
		prot1 = sps[pair[0]] + pair[2][0]
		prot2 = sps[pair[1]] + pair[2][1]
		if prot1 not in prot_ind and prot2 not in prot_ind:
			clusters.append([])
			clusters[cnt].append(prot1)
			clusters[cnt].append(prot2)
			prot_ind[prot1] = cnt
			prot_ind[prot2] = cnt
			cnt+=1
		else:
			if prot1 in prot_ind and prot2 not in prot_ind:
				ind1 = prot_ind[prot1]
				if len(clusters[ind1]) < 6:
					prot_ind[prot2] = ind1
					clusters[ind1].append(prot2)
			if prot1 not in prot_ind and prot2 in prot_ind:
				ind2 = prot_ind[prot2]
				if len(clusters[ind2]) < 6:
					prot_ind[prot1] = ind2
					clusters[ind2].append(prot1)"""

	return new_clusters


def seed_gen(k, filename, pairs_dict, letters):
	all_sp = []
	sps = []
	sps_int = dict()
	for i in range(k):
		all_sp.append(i)
		sps.append(letters[i])
		sps_int[letters[i]] = i
	pairs = []
	for i in range(k):
		for j in range(i+1, k):
			for pair in pairs_dict[i][j]:
				if pairs_dict[i][j][pair] > 40:
					pairs.append([i, j, pair, pairs_dict[i][j][pair]])
	pairs.sort(key = lambda tup: - 1 * tup[3])
	matched_pairs = []
	matched_sps = []
	vec_prot = []
	for i in range(k):
		matched_pairs.append(defaultdict(list))
		matched_sps.append(defaultdict(set))
		vec_prot.append('-1')
	all_matched = []
	cnt = 0
	matched_prot = dict()
	for run in range(1):
		for pair in pairs:
			sp1 = pair[0]
			sp2 = pair[1]
			prot1 = pair[2][0]
			prot2 = pair[2][1]

			sp_prot1 = sps[sp1] + prot1
			sp_prot2 = sps[sp2] + prot2
			if sp_prot1 in matched_prot:
				ind1 = matched_prot[sp_prot1]
				if sp_prot2 not in matched_prot:
					if all_matched[ind1][sp2] == '-1':
						all_matched[ind1][sp2] = prot2
						matched_prot[sp_prot2] = ind1
				else:
					ind2 = matched_prot[sp_prot2]
					######################################################################################
					"""if all_matched[ind1][sp2] == '-1' and all_matched[ind2][sp1] == '-1':
						merge_flag = True
						cnt = 0
						for index in range(k):
							if all_matched[ind1][index] != '-1' and all_matched[ind2][index] != '-1':
								cnt+=1
						if cnt == -1:
							merge_flag = False
						#merge_flag = True
						if ind1 < ind2:
							tind = ind1
							ind1 = ind2
							ind2 = tind						
						if merge_flag and ind1 != ind2:
							#print all_matched[ind2], all_matched[ind1], cnt
							for index in range(k):
								if all_matched[ind2][index] != '-1' and all_matched[ind1][index] == '-1':
									all_matched[ind1][index] = all_matched[ind2][index]
									sp_tmp = sps[index] + all_matched[ind2][index]
									matched_prot[sp_tmp] = ind1
									all_matched[ind2][index] = '-1'
							#print all_matched[ind2], all_matched[ind1]
							if cnt == -1:
								print all_matched[ind2], all_matched[ind1]
								for index in range(k):
									if all_matched[ind2][index] != '-1':
										sp_tmp = sps[index] + all_matched[ind2][index]
										del matched_prot[sp_tmp]
										all_matched[ind2][index] = '-1'
								print all_matched[ind2], all_matched[ind1]
								print "+++++++++++++++++++++++++++++++++++"""
					######################################################################################
					"""if all_matched[ind1][sp2] == '-1' and all_matched[ind2][sp1] == '-1':
						if ind1 < ind2:
							tind = ind1
							ind1 = ind2
							ind2 = tind
						if ind1 != ind2:
							for index in range(k):
								if all_matched[ind1][index] == '-1':
									if all_matched[ind2][index] != '-1':
										all_matched[ind1][index] = all_matched[ind2][index]
										sp_tmp = sps[index] + all_matched[ind2][index]
										matched_prot[sp_tmp] = ind1
										all_matched[ind2][index] = '-1'
							cnt = 0
							ind0 = -1
							for index in range(k):
								if all_matched[ind2][index] != '-1':
									cnt+=1
									ind0 = index
							if cnt == 1:
								node = all_matched[ind2][ind0]
								sp_tmp = sps[ind0] + all_matched[ind2][ind0]
								all_matched[ind2][ind0] = '-1'
								del matched_prot[sp_tmp]
								#print all_matched[ind2], all_matched[ind1], cnt, all_matched[ind2][ind0], sp_tmp, matched_prot[sp_tmp] """
					if all_matched[ind1][sp2] == '-1' and all_matched[ind2][sp1] == '-1':
						if ind1 != ind2:
							if ind1 > ind2:
								for index in range(k):
									if all_matched[ind1][index] == '-1':
										if all_matched[ind2][index] != '-1':
											all_matched[ind1][index] = all_matched[ind2][index]
											sp_tmp = sps[index] + all_matched[ind2][index]
											matched_prot[sp_tmp] = ind1
											all_matched[ind2][index] = '-1'
							else:
								for index in range(k):
									if all_matched[ind2][index] == '-1':
										if all_matched[ind1][index] != '-1':
											all_matched[ind2][index] = all_matched[ind1][index]
											sp_tmp = sps[index] + all_matched[ind1][index]
											matched_prot[sp_tmp] = ind2
											all_matched[ind1][index] = '-1'
							#print all_matched[ind2], all_matched[ind1]
			else:
				if sp_prot2 in matched_prot:
					ind2 = matched_prot[sp_prot2]
					if all_matched[ind2][sp1] == '-1':
						all_matched[ind2][sp1] = prot1
						matched_prot[sp_prot1] = ind2
				else:
					all_matched.append(copy.copy(vec_prot))
					ind = len(all_matched) - 1
					all_matched[ind][sp1] = prot1
					all_matched[ind][sp2] = prot2
					matched_prot[sp_prot1] = ind
					matched_prot[sp_prot2] = ind
	tot = 0

	qunits_size = []
	quints = []
	all_prots = set()
	for vals in all_matched:
		cnt = 0
		for i in range(k):
			if vals[i] != '-1':
				cnt +=1
		if cnt > 1:
			for i in range(k):
				if vals[i] != '-1':
					prot_sp = sps[i] + vals[i]
					all_prots.add(prot_sp)
			tot +=1
			quints.append(vals)
			qunits_size.append([vals, cnt])


	sets = []
	for i in range(int(k)):
		sets.append(set())
	fto = open(filename, 'w')
	#print len(quints)
	for quint in reversed(quints):
		prots = []
		for prot in quint:
			prots.append(prot)
		flag = True
		for i in range(int(k)):
			if prots[i] in sets[i] and prots[i] != '-1':
				flag = False
		if flag:
			for i in range(int(k)):
				if prots[i] != '-1':
					sets[i].add(prots[i])

			for i in range(int(k)):
				fto.write(prots[i])
				if i < int(k) - 1:
					fto.write('\t')
				else:
					fto.write('\n')
	fto.close()
	"""for i in range(int(k)):
		print len(sets[i]),
	print """
	return quints

def load_synth(k, alg, family, min_blast):
	families = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
	all_graphs = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
	all_letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
	graphs = all_graphs[:int(k)]

	directory='../NAPAbench/'+ str(k) + '-way/' + alg + '_set/Family_' + family + '/'
	
	dict_int = dict()
	for i in range(int(k)):
		dict_int[graphs[i]] = i
	letters = all_letters[:int(k)]

	

	vector = []
	for i in range(int(k)):
		vector.append([])
	lr = []
	rl = []
	pairs_dict = []
	nodes = []
	for i in range(int(k)):
		lr.append(copy.copy(vector))
		rl.append(copy.copy(vector))
		pairs_dict.append(copy.copy(vector))
	for i in range(int(k)):
		for j in range(int(k)):
			lr[i][j] = defaultdict(list)
			rl[i][j] = defaultdict(list)
			pairs_dict[i][j] = dict()
	for i in range(int(k)):
		for j in range(i+1, int(k)):
			filename = directory + graphs[i] + '-' + graphs[j]+ '.sim'
			f = open(filename, 'r')
			for line in f:
				#print parsed
				parsed = line.replace('\n', '').split()
				#print parsed[0], parsed[1], filename
				node1 = parsed[0]
				node2 = parsed[1]
				sim = parsed[2]
				if float(sim) >= min_blast:
					pairs_dict[i][j][(node1, node2)] = float(sim)
					pairs_dict[j][i][(node2, node1)] = float(sim)
					lr[i][j][node1].append((node2, float(sim)))
					rl[j][i][node2].append((node1, float(sim)))

	return pairs_dict, lr, rl


def is_number(s):
    try:
        return True, float(s)
    except ValueError:
        return False, -1