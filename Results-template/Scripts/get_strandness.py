import sys

# This script was written specifically for RNASeqTechDev
# returns strandness
# 0 --> unstranded
# 1 --> forward stranded
# 2 --> reverse stranded

def most_common(lst):
    return max(set(lst), key=lst.count)

def get_correct_strandness(basefolder,g):
    glines=list(map(lambda x:x.strip(),open(basefolder+"/"+g+".RnaSeqMetrics.txt")))
    for ln,l in enumerate(glines):
        if l.startswith("PF_BASES"):
            headers=l.split("\t")
            break
    values=glines[ln+1].split()
    strandness=float(values[headers.index("PCT_CORRECT_STRAND_READS")])
    if strandness < 0.25:
        return 1 # sense
    elif strandness <= 0.75 :
        return 0 # unstranded
    else:
        return 2 # antisense

basefolder=sys.argv[1] # folder that contains groups.tab file
groupsfile=basefolder+"/groups.tab"
groups=filter(None,map(lambda x:x.strip().split(),open(groupsfile).readlines()))
strandnesses=[]
for g in groups:
    strandnesses.append(get_correct_strandness(basefolder,g[0]))

print(most_common(strandnesses))
