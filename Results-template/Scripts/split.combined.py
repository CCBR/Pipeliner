#! /usr/bin/env python

F=open("combined.gvcf","r")
snp=open("all.snp.vcf","w")
indel=open("all.indel.vcf","w")

for line in F.read().split("\n"):
    if len(line)==0:
        break
    if line.startswith("#"):
        snp.write(line+"\n")
        indel.write(line+"\n")

    else:
        L=line.split("\t")
        if len(L[3])==1 & len(L[4])==1:
            snp.write(line+"\n")
        else:
            indel.write(line+"\n")

F.close()
snp.close()
indel.close()

