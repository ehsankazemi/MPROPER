import numpy as np


filename = "scripts/150.txt"
filename = "out.txt"
sizes = np.zeros(6)
f = open(filename, 'r')
for line in f:
	parsed = line.strip().split()
	cnt = 0
	for val in parsed:
		if val != "-1":
			cnt+=1
	sizes[cnt] +=1

print(sizes)