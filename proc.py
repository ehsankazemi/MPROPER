import numpy as np

sps = ['ce', 'dm', 'hs', 'mm', 'sc']


for i in range(5):
	for j in range(i + 1, 5):
		sp1 = sps[i]
		sp2 = sps[j]
		filenameR = 'sorted_blast_' + sp1 + '-' + sp2 + '.blast'
		filenameW =  'MPROPER/dataset/sequence_similarity/' +sp1 + '-' + sp2 + '.blast'
		f = open(filenameR, 'r')
		fto = open(filenameW, 'w')
		for line in f:
			parsed = line.strip().split()
			if int(parsed[6]) >= 40:
				fto.write(parsed[0] + '\t' + parsed[3] + '\t' + parsed[6] + '\n')
			else:
				break
		fto.close()
