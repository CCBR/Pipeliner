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

filterType = ['nFeature','nCount','mitoPct']
filtStatus = ['preFilter','postFilter']

filtList = []
for q in filterType:
	for r in filtStatus:
		filtList.append(str(q+"_"+r))

contList = []
contrasts = open("contrasts.tab")
for line in contrasts:
	line = line.strip("\n")
#	print(re.sub('\t', '-',  line))
	contList.append(re.sub('\t', '-',  line))

print("{}".format(contList[0]))

resList = resolution.split(",")

intType = ["seurat_batch_corrected","merged"]

rule all:
	params:
		batch='--time=168:00:00',
	input:
		# expand(join(workpath,"QC","{name}.RData"),name=samples),
		#QC file outputs
		expand(join(workpath,"filtered","{name}.rds"),name=samples),
		expand(join(workpath,"flags","{name}.txt"),name=samples),
		expand(join(workpath,"filtered","{name}_doublets.rds"),name=samples),
		#QC image outputs
		expand(join(workpath,"QC","{name}","images","filterStats_{name}.png"),name=samples),
		expand(join(workpath,"QC","{name}","images","cellsRemovedVenn_{name}.png"),name=samples),
		expand(join(workpath,"QC","{name}","images","{filt}_{time}_{name}.pdf"),name=samples,filt=filterType,time=filtStatus),
		expand(join(workpath,"QC","{name}","images","clusterResolution_{res}_{name}.pdf"),name=samples,res=resList),
		expand(join(workpath,"QC","{name}","images","silhouetteResolution_{res}_{name}.pdf"),name=samples,res=resList),
		expand(join(workpath,"QC","{name}","images","primaryAnnotation_{name}.pdf"),name=samples),
		expand(join(workpath,"QC","{name}","images","doublets_{name}.pdf"),name=samples),
		expand(join(workpath,"QC","{name}","images","cellCycle_{name}.pdf"),name=samples),

		#join(workpath,"QC","samples_QC.html"),
		expand(join(workpath,"QC","QC_Report_{name}.html"),name=samples),
		expand(join(workpath,"integration","seurat_batch_corrected","{myGroups}","{myGroups}.rds"),myGroups=contList),
		expand(join(workpath,"integration","merged","{myGroups}","{myGroups}.rds"),myGroups=contList),
		expand(join(workpath,"QC","integration","{myGroups}"),myGroups=contList),
		
rule qc_scrna:
	input: 
		join(workpath,"{name}.h5"),
	output: 
		rds=join(workpath,"filtered","{name}.rds"),
		#rdata=join(workpath,"QC","{name}.RData"),
		qc = join(workpath,"flags","{name}.txt"),
		doublets=join(workpath,"filtered","{name}_doublets.rds"),
		filtStats=join(workpath,"QC","{name}","images","filterStats_{name}.png"),
		venn=join(workpath,"QC","{name}","images","cellsRemovedVenn_{name}.png"),
		filtVln=expand(join(workpath,"QC","{{name}}","images","{filt}_{{name}}.pdf"),filt=filtList),
		clustUMAP=expand(join(workpath,"QC","{{name}}","images","clusterResolution_{res}_{{name}}.pdf"),res=resList),
		clustSil=expand(join(workpath,"QC","{{name}}","images","silhouetteResolution_{res}_{{name}}.pdf"),res=resList),
		annotUMAP=join(workpath,"QC","{name}","images","primaryAnnotation_{name}.pdf"),
		doubletUMAP=join(workpath,"QC","{name}","images","doublets_{name}.pdf"),
		cellCycleUMAP=join(workpath,"QC","{name}","images","cellCycle_{name}.pdf"),
		imageDir = join(workpath,"QC","{name}","images"),
	params:
		rname='pl:qc_scrna',
		specie = specie, 
		resolution = resolution,
		clustAlg = clustAlg,
		annotDB = annotDB,
		citeseq = citeseq,
		outDir= join(workpath,"QC"),
		imageDir = join(workpath,"QC","{name}","images")

	shell: """

module load R/3.6.1;
Rscript Scripts/scrnaQC.R {input} {output.rds} {params.imageDir} {params.specie} {params.resolution} {params.clustAlg} {params.annotDB} {params.citeseq};
touch {output.qc}
        """
#Rscript Scripts/scrnaQC.R {input} {output.rds} {output.rdata} {params.specie} {params.resolution} {params.clustAlg} {params.annotDB} {params.citeseq};

rule qcReport_scrna:
	input:
		# files = expand(join(workpath,"QC","{name}.RData"),name=samples),
		path=join(workpath,"QC","{name}","images"),
	output:
		out=join(workpath,"QC","QC_Report_{name}.html"),
	params:
		rname='pl:qcReport_scrna',
		resolution=resolution,
		sample = "{name}",
		outPath = join(workpath,"QC")
	shell: """
module load R/3.6.1;
Rscript Scripts/scrnaQC_Reports.R {input.path} {params.outPath} {params.resolution} {params.sample};
	"""
# mv Scripts/samples_QC.html QC
	
rule integratedBatch:
	input:
		rds=expand(join(workpath,"filtered","{name}.rds"),name=samples),
		flag=expand(join(workpath,"flags","{name}.txt"),name=samples)
	output:
		rdsBatch=join(workpath,"integration","seurat_batch_corrected","{myGroups}","{myGroups}.rds"),
		mergeRDS=join(workpath,"integration","merged","{myGroups}","{myGroups}.rds"),
		jointImageDir=join(workpath,"QC","integration","{myGroups}","images")
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
Rscript Scripts/integrateBatches.R {params.dir} {output.rdsBatch} {output.mergeRDS} {output.jointImageDir} {params.specie} {params.resolution} {params.clustAlg} {params.annotDB} {params.nAnchors} {params.citeseq} {params.groups} {params.contrasts};
	"""



