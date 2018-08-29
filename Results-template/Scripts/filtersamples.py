import sys,os
folderpath=sys.argv[1]
rver=sys.argv[2]
rscript=sys.argv[3]
sampletable=sys.argv[4]
rawcountsfile=sys.argv[5]
newsampletable=sys.argv[6]
filteredcountstable=sys.argv[7]

degdir=folderpath.strip().split("/")[-1]

v1=degdir.split("DEG_")[1]
minsamples=v1.split("_")[-1]
cpm_cutoff=v1.split("_")[-2]

v2=v1.split("_")
v2.pop(-1)
v2.pop(-1)
v3="_".join(v2)
g1=v3.split("-")[0]
g2=v3.split("-")[1]

print(" ".join(["module load",rver,"&&","Rscript",rscript,folderpath,g1,g2,cpm_cutoff,minsamples,sampletable,rawcountsfile,newsampletable,filteredcountstable]))

os.system(" ".join(["module load",rver,"&&","Rscript",rscript,folderpath,g1,g2,cpm_cutoff,minsamples,sampletable,rawcountsfile,newsampletable,filteredcountstable]))
