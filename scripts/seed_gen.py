import numpy as np
import copy
from collections import defaultdict
import json
import pickle
import sys

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
	return quints

def read_blast(sp1, sp2, min_blast):
	f = open( '../dataset/sequence_similarity/' + sp1 +'-' + sp2 +'.blast', 'r' )
	lr = defaultdict(list)
	rl = defaultdict(list)
	pairs_dict12 = dict()
	pairs_dict21 = dict()
	for line in f:
		parsed = line.replace( "\n", "" ).split()
		if int(parsed[2]) >= min_blast:
			prot1 = parsed[0]
			prot2 = parsed[1]
			pairs_dict12[(prot1, prot2)] = int(parsed[2])
			pairs_dict21[(prot2, prot1)] = int(parsed[2])
			lr[prot1].append((prot2, int(parsed[2])))
			rl[prot2].append((prot1, int(parsed[2])))	
		else:
			break
	print len(lr), min_blast
	return lr, rl, pairs_dict12, pairs_dict21

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


filename = sys.argv[1]
min_blast = int(sys.argv[2])
k = 5 #number of species
sps = ['ce', 'dm', 'hs', 'mm', 'sc'] #list of species
lr, rl, pairs_dict = load_all_data(min_blast, sps) #read the blast bit-score similarities
seeds = seed_gen(int(k), filename, pairs_dict, sps) #generate the set of seed tuples

