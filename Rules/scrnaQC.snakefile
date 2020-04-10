from snakemake.utils import R
from os.path import join
from os import listdir
import os
import re

configfile: "run.json"
samples=config['project']['groups']['rsamps']
workpath=config['project']['workpath']

# workpath="/data/abdelmaksoudaa/klion"
# samples = [f for f in os.listdir('.') if re.search('4H', f)]

specie = config['project']['SPECIES'] #human
resolution = config['project']['CLUSTRESOLUTION']#"0.8"
clustAlg = config['project']['CLUSTALG'] #1-3
annotDB = config['project']['ANNOTDB']#"HPCA" #mouse
citeseq = config['project']['CITESEQ']
nAnchors = "2000"
groups = "YES" #YES/NO


contList = []
contrasts = open("contrasts.tab")
for line in contrasts:
	line = line.strip("\n")
#	print(re.sub('\t', '-',  line))
	contList.append(re.sub('\t', '-',  line))


intType = ["seurat_batch_corrected","merged"]

rule all:
	params:
		batch='--time=168:00:00',
	input:
		expand(join(workpath,"QC","{name}.RData"),name=samples),
		expand(join(workpath,"filtered","{name}.rds"),name=samples),
		expand(join(workpath,"flags","{name}.txt"),name=samples),
		join(workpath,"QC","samples_QC.html"),		
		expand(join(workpath,"integration","seurat_batch_corrected","{myGroups}","{myGroups}.rds"),myGroups=contList),
		expand(join(workpath,"integration","merged","{myGroups}","{myGroups}.rds"),myGroups=contList),
		
rule qc_scrna:
	input: 
		join(workpath,"{name}.h5")
	output: 
		rds=join(workpath,"filtered","{name}.rds"),
		rdata=join(workpath,"QC","{name}.RData"),
		qc = join(workpath,"flags","{name}.txt"),
	params:
		rname='pl:qc_scrna',
		specie = specie, 
		resolution = resolution,
		clustAlg = clustAlg,
		annotDB = annotDB,
		citeseq = citeseq,
		outDir= join(workpath,"QC")
	shell: """

module load R/3.6.1;
Rscript Scripts/scrnaQC.R {input} {output.rds} {output.rdata} {params.specie} {params.resolution} {params.clustAlg} {params.annotDB} {params.citeseq};
touch {output.qc}
        """

rule qcReport_scrna:
	input:
		files = expand(join(workpath,"QC","{name}.RData"),name=samples),
		path=join (workpath, "QC")
	output:
		join(workpath,"QC","samples_QC.html")
	params:
		rname='pl:qcReport_scrna',
	shell: """
module load R/3.6.1;
Rscript Scripts/scrnaQC_Reports.R {input.path};
mv Scripts/samples_QC.html QC
	"""
	
rule integratedBatch:
	input:
		rds=expand(join(workpath,"filtered","{name}.rds"),name=samples),
		flag=expand(join(workpath,"flags","{name}.txt"),name=samples)
	output:
		rdsBatch=join(workpath,"integration","seurat_batch_corrected","{myGroups}","{myGroups}.rds"),
		mergeRDS=join(workpath,"integration","merged","{myGroups}","{myGroups}.rds"),
	params:
		rname='pl:integratedBatch',
		batch='--cpus-per-task=8 --mem=48g --time=24:00:00',
		specie = specie,
		resolution = resolution,
		clustAlg = clustAlg,
		annotDB = annotDB,
#		outDirSeurat= rdsBatch#join(workpath,"integration","seurat_batch_corrected_RDS"),
#		outDirMerged= mergeRDS#join(workpath,"integration","mergedRDS"),
		nAnchors = nAnchors,
		groups = groups,
		dir = join(workpath,"filtered"),
		citeseq = citeseq,
		contrasts = "{myGroups}"
	shell: """
module load R/3.6.1;
Rscript Scripts/integrateBatches.R {params.dir} {output.rdsBatch} {output.mergeRDS} {params.specie} {params.resolution} {params.clustAlg} {params.annotDB} {params.nAnchors} {params.citeseq} {params.groups} {params.contrasts};
	"""



