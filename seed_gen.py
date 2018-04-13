import numpy as np
import copy
from collections import defaultdict
import json
import pickle
import utils
name= 'EAL5'
gos = pickle.load(open("GO"+ name +".dat", "rb"))
k_tmp='5'
filename = 'test.txt'
min_blast = 100
k = 5
sps = ['ce', 'dm', 'hs', 'mm', 'sc']
lr, rl, pairs_dict = utils.load_all_data(min_blast, sps)
letters = sps
filename = 'test.txt'
k = int(k)
annotated = 0
correct = 0
c_clustres = np.zeros(int(k) + 1)
nb_c = 0
seeds = utils.seed_gen(int(k), filename, pairs_dict, letters)

exit()
for seed in seeds:
	nodes = []
	for i in range(k):
		if seed[i] != '-1':
			nodes.append(letters[i] + seed[i])
	cnt = 0
	for node in nodes:
		if node in gos:
			go = gos[node]
			if cnt == 0:
				inter_sec = go
			else:
				inter_sec = utils.intersect(inter_sec, gos[node])
			cnt += 1
	if cnt == len(seed) or cnt > 1:
		annotated += 1
		if len(inter_sec) > 0:
			correct += 1.0
			nb_sp = set()
			for node in seed:
				nb_sp.add(node[0])
			c_clustres[len(nb_sp)]+=1
			nb_c += len(seed)
print annotated, correct, correct / annotated, c_clustres, nb_c
