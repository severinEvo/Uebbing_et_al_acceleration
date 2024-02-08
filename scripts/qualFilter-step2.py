#!/usr/bin/env python
import gzip
import sys

use = sys.argv[1] # "qualFilter_useTips.bed.gz"
clades = sys.argv[2] # "annotations/clalist.pyData"
outfile = sys.argv[3] # "qualFilter_useClades.tsv"

useDict = {}
with gzip.open(use,'r') as fin:
	for line in fin:
		line = bytes.decode(line).rstrip('\r\n')
		s = line.split("\t")
		pos = s[0] + ":" + s[1] + "-" + s[2]
		if pos not in useDict:
			useDict[pos] = [s[3]]
		else:
			useDict[pos].append(s[3])

exec(open(clades).read())
fout = open(outfile,"w")
for clade in clalist:
	in1_list = clalist[clade]['in1_list']
	in2_list = clalist[clade]['in2_list']
	for p in useDict:
		resIn1 = any(item in useDict[p] for item in in1_list)
		resIn2 = any(item in useDict[p] for item in in2_list)
		if resIn1 and resIn2:
			fout.write(p + "\t" + clade + "\n")

fout.close()
