from __future__ import print_function
import os

files = list(filter(lambda x:"HPC_usage_table.txt" in x,os.listdir("./")))
#print(files)

lines = list()
for f in files:
	#ignore existing "._" files created by sublime text 
	if not f.startswith('._'):
		lines.extend(open(f).readlines())
outlines=list()
outlines.append(lines.pop(0))
lines=list(filter(lambda x:not x.startswith("Job"),lines))
lines_set=set(lines)
for l in lines:
	if lines.count(l)==1:
		outlines.append(l)
for f in outlines:
	print(f.strip())
