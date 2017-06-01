import sys,re,pprint
import numpy as np
pp = pprint.PrettyPrinter(indent=4)

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


deseqFCfile=sys.argv[1]
idpos = int(sys.argv[2])
logfcpos = int(sys.argv[3])
species = sys.argv[4]

logfcdata={re.sub("\"","",a[idpos]):float(a[logfcpos]) for a in list(map(lambda x:x.strip().split("\t"),open(deseqFCfile).readlines()[1:]))}

#pp.pprint(logfcdata)

if "hg" in species:
	chr_numbers=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

if "mm" in species:
	chr_numbers=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', 'X', 'Y']

chrs=[ "chr"+a for a in chr_numbers ]

binsize=50
list_genes=1 #Changeable
for chrom in chrs:
	bedfile="../bedfiles/karyobed."+chrom+".bed"
	lines=list(map(lambda x:x.strip().split("\t"),open(bedfile).readlines()))
	genelist={a[4]:"###".join([ a[0], a[1], a[2]]) for a in lines}
	datapoints=[]
	for line in lines:
		gene=line[4]
		if gene in logfcdata:
			logfc=logfcdata[gene]
			#if logfc<-3:
			#	logfc=-3
			#if logfc>3:
			#	logfc=3
			datapoints.append([chrom,line[1],line[2],logfc,gene])
	#multiple=int(len(datapoints)/binsize)
	#print(len(datapoints))
	#datapoints=datapoints[:multiple*binsize]
	#print(len(datapoints))
	for databin in list(chunks(datapoints,binsize)):
		fc=str(np.mean(np.array([a[3] for a in databin])))
		start=databin[0][1]
		end=databin[-1][2]
		if list_genes==1:
			genes=",".join([a[4] for a in databin])
		else:
			genes="."
		print("\t".join([chrom,start,end,genes,fc]))
		#print("\t".join([chrom,start,end,".",fc]))
	#fcavg=np.mean(fc.reshape(-1, binsize), axis=1)
	#print(fcavg)
	#pp.pprint(genelist)
