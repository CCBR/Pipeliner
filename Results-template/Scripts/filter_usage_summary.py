import os
files = list(filter(lambda x:"HPC_usage_table.txt" in x,os.listdir("./")))
lines = list()
for f in files:
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
