import pandas as pd
import os,sys
from os.path import join
import numpy

sampletable=pd.read_csv(sys.argv[3],header=0,sep="\t")
sampletable['sampleName']=sampletable['sampleName'].map(lambda x:x.split("/")[-1].split(".star.count.")[0])

outdir=sys.argv[1]
igfiles=sys.argv[2].split()
contrasts=sys.argv[4].split()
rsemref=sys.argv[5]
itype=sys.argv[6]

annotate=pd.read_csv(sys.argv[7],header=None,sep=" ",usecols=[0,2])
annotate.columns=["transcriptID","geneID"]

condition2sampleName=dict()
sampleName2igfile=dict()
for i,sample in sampletable.iterrows():
    if not sample['condition'] in condition2sampleName:
        condition2sampleName[sample['condition']]=list()
    condition2sampleName[sample['condition']].append(sample['sampleName'])
    sampleName2igfile[sample['sampleName']]=list(filter(lambda x:sample['sampleName']+".RSEM" in x,igfiles))[0]

# print(sampletable)
# print(condition2sampleName)
# print(sampleName2igfile)

for i in range(0,len(contrasts),2):
    pair=contrasts[i]+"_vs_"+contrasts[i+1]
    cmfile=join(outdir,".".join([pair,itype,"counts","matrix"]))
    ebseqfile=join(outdir,".".join([pair,itype,"EBSeq"]))
    files=list()
    files.extend(list(map(lambda x:sampleName2igfile[x],condition2sampleName[contrasts[i]])))
    files.extend(list(map(lambda x:sampleName2igfile[x],condition2sampleName[contrasts[i+1]])))
    cmd=list()
    cmd.append("rsem-generate-data-matrix")
    cmd.extend(files)
    cmd.append(">")
    cmd.append(cmfile)
    cmd=" ".join(cmd)
    os.system(cmd)
    print('Running CMD: {}'.format(cmd))
    cmd=list()
    cmd.append("rsem-run-ebseq")
    if itype=="isoform":
    	cmd.append("--ngvector")
    	cmd.append(rsemref+".transcripts.ngvec")
    cmd.append(cmfile)
    cmd.append("%d,%d"%(len(condition2sampleName[contrasts[i]]),len(condition2sampleName[contrasts[i+1]])))
    cmd.append(ebseqfile)
    cmd=" ".join(cmd)
    os.system(cmd)
    print('Running CMD: {}'.format(cmd))
    x=pd.read_csv(ebseqfile,sep="\t",header=0)
    x.index.name="transcriptID"
    x.reset_index(inplace=True)
    x=pd.merge(annotate,x,on="transcriptID")
    x.sort_values("PPDE",inplace=True,ascending=False)
    x["logPostFC"]=x["PostFC"].apply(numpy.log2)
    x.to_csv(ebseqfile,sep="\t",index=False)

finaloutfile=open(join(outdir,"EBSeq_%s_completed.txt"%(itype)),'w')
finaloutfile.write("%s"%str(sampletable))
finaloutfile.close()
