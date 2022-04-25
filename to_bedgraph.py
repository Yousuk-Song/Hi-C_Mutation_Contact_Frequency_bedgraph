#!/home/eaststar0/miniconda3/bin/python

import sys
dfh=open(sys.argv[1],'r')
rfh=open(sys.argv[1][:-4]+'.bedgraph','w')
window = int(sys.argv[2])

count=0
D = {}
for line in dfh:
	line = line.rstrip()
	if line not in D:
		D[line] = 1
	else:
		D[line] += 1

for key in D:
	[chrom, pos] = key.split()
	rfh.write(f"{chrom}\t{pos}\t{int(pos) + window}\t{D[key]}\n")

rfh.close()

