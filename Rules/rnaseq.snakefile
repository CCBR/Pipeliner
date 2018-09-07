from snakemake.utils import R
from os.path import join
from os import listdir
import os

def check_existence(filename):
  if not os.path.exists(filename):
    exit("File: %s does not exists!"%(filename))

def check_readaccess(filename):
  check_existence(filename)
  if not os.access(filename,os.R_OK):
    exit("File: %s exists, but cannot be read!"%(filename))

def check_writeaccess(filename):
  check_existence(filename)
  if not os.access(filename,os.W_OK):
    exit("File: %s exists, but cannot be read!"%(filename))

def createConstrasts(cList):
  contrastsList = []
  for i in range(0, len(cList)-1, 2):
    contrastsList.append("-".join(cList[i:i+2]))
  return contrastsList

def get_cpm_cutoff(degdir,contrasts_w_cpm_cutoff_list,cpm_cutoff_list):
	return cpm_cutoff_list[contrasts_w_cpm_cutoff_list.index(degdir)]

configfile: "run.json"

samples=config['project']['groups']['rsamps']

workpath=config['project']['workpath']

contrastsList = createConstrasts(config['project']['contrasts']['rcontrasts'])
cpm_cutoff_list = config['project']['contrasts']['rcontrasts_cpm_cutoff']
mincounts = config['project']['contrasts']['rcontrasts_min_counts']

contrasts_w_cpm_cutoff_list = []
for i in zip(contrastsList,cpm_cutoff_list,mincounts):
	contrasts_w_cpm_cutoff_list.append(str(i[0])+"_"+str(i[1])+"_"+str(i[2]))

se=""
pe=""
workpath = config['project']['workpath']

if config['project']['nends'] == 2 :
  pe="yes"
elif config['project']['nends'] == 1 :
  se="yes"

check_readaccess(join(workpath,"groups.tab"))
check_readaccess(join(workpath,"contrasts.tab"))

trim_dir='trim'
star_dir="STAR_files"
bams_dir="bams"
log_dir="logfiles"
rseqc_dir="RSeQC"
# rsemi_dir="DEG_RSEM_isoforms"
# rsemg_dir="DEG_RSEM_genes"
# salmon_dir="Salmon"
# subreadg_dir="DEG_Subread_genes"
# subreadj_dir="DEG_Subread_junctions"
# subreadgj_dir="DEG_Subread_genejunctions"



if config['project']['DEG'] == "yes" and config['project']['TRIM'] == "yes":
	rule all:
		params: 
			batch="--time=168:00:00",
		
		
		input: 
			expand(join(workpath,"DEG_{deg_dir}","EBSeq_isoform_completed.txt"),deg_dir=contrasts_w_cpm_cutoff_list),
			expand(join(workpath,"DEG_{deg_dir}","EBSeq_gene_completed.txt"),deg_dir=contrasts_w_cpm_cutoff_list),
			expand(join(workpath,"DEG_{deg_dir}","RawCountFile_RSEM_genes_filtered.txt"), deg_dir=contrasts_w_cpm_cutoff_list),
			expand(join(workpath,"DEG_{deg_dir}","DESeq2_PCA.png"),deg_dir=contrasts_w_cpm_cutoff_list),
			expand(join(workpath,"DEG_{deg_dir}","edgeR_prcomp.png"),deg_dir=contrasts_w_cpm_cutoff_list),
			expand(join(workpath,"DEG_{deg_dir}","limma_MDS.png"),deg_dir=contrasts_w_cpm_cutoff_list),
# 			expand(join(workpath,"DEG_{deg_dir}","PcaReport.html"),deg_dir=contrasts_w_cpm_cutoff_list),
# 			expand(join(workpath,rsemg_dir,"limma_edgeR_DESeq2_{con}_vennDiagram.png"),con=contrastsList),


rule deseq2:
	input:
		file1=join(workpath,"DEG_{degdir}","sampletable.txt"),
		file2=join(workpath,"DEG_{degdir}","RawCountFile_RSEM_genes_filtered.txt"),
	output:
		join(workpath,"DEG_{degdir}","DESeq2_PCA.png"),
	params:
		rname='pl:deseq2',
		batch='--mem=24g --time=10:00:00',
		outdir=join(workpath,"DEG_{degdir}"),
		contrasts=" ".join(config['project']['contrasts']['rcontrasts']),
		dtype="RSEM_genes",
		refer=config['project']['annotation'],
		projectId=config['project']['id'],
		projDesc=config['project']['description'],
		gtffile=config['references'][pfamily]['GTFFILE'],
		karyobeds=config['references'][pfamily]['KARYOBEDS'],
		rver=config['bin'][pfamily]['tool_versions']['RVER'],
		rscript1=join(workpath,"Scripts","deseq2call.R"),
		rscript2=join(workpath,"Scripts","Deseq2Report.Rmd"),
	shell: """
cp {params.rscript1} {params.outdir}
cp {params.rscript2} {params.outdir}
cd {params.outdir}
module load {params.rver}
Rscript deseq2call.R '{params.outdir}' '{input.file1}' '{input.file2}' '{params.contrasts}' '{params.refer}' '{params.projectId}' '{params.projDesc}' '{params.gtffile}' '{params.dtype}' '{params.karyobeds}'
"""


rule edgeR:
  input:
    file1=join(workpath,"DEG_{degdir}","sampletable.txt"),
    file2=join(workpath,"DEG_{degdir}","RawCountFile_RSEM_genes_filtered.txt"),
  output:
    join(workpath,"DEG_{degdir}","edgeR_prcomp.png"),
  params: 
    rname='pl:edgeR',
    batch='--mem=24g --time=10:00:00',
    outdir=join(workpath,"DEG_{degdir}"),
    contrasts=" ".join(config['project']['contrasts']['rcontrasts']),
    dtype="RSEM_genes",
    refer=config['project']['annotation'],
    projectId=config['project']['id'],
    projDesc=config['project']['description'],
    gtffile=config['references'][pfamily]['GTFFILE'],
    karyobeds=config['references'][pfamily]['KARYOBEDS'],
    rver=config['bin'][pfamily]['tool_versions']['RVER'],
    rscript1=join(workpath,"Scripts","edgeRcall.R"),
    rscript2=join(workpath,"Scripts","EdgerReport.Rmd"),
  shell: """
cp {params.rscript1} {params.outdir}
cp {params.rscript2} {params.outdir}
cd {params.outdir}
module load {params.rver}
Rscript edgeRcall.R '{params.outdir}' '{input.file1}' '{input.file2}' '{params.contrasts}' '{params.refer}' '{params.projectId}' '{params.projDesc}' '{params.gtffile}' '{params.dtype}' '{params.karyobeds}'
"""

rule limmavoom:
  input: 
    file1=join(workpath,"DEG_{degdir}","sampletable.txt"),
    file2=join(workpath,"DEG_{degdir}","RawCountFile_RSEM_genes_filtered.txt"),
  output: 
    join(workpath,"DEG_{degdir}","limma_MDS.png")
  params:
    rname='pl:limmavoom',
    batch='--mem=24g --time=10:00:00',
    outdir=join(workpath,"DEG_{degdir}"),
    contrasts=" ".join(config['project']['contrasts']['rcontrasts']),
    dtype="RSEM_genes",
    refer=config['project']['annotation'],
    projectId=config['project']['id'],
    projDesc=config['project']['description'],
    gtffile=config['references'][pfamily]['GTFFILE'],
    karyobeds=config['references'][pfamily]['KARYOBEDS'],
    rver=config['bin'][pfamily]['tool_versions']['RVER'],
    rscript1=join(workpath,"Scripts","limmacall.R"),
    rscript2=join(workpath,"Scripts","LimmaReport.Rmd"),
  shell: """
cp {params.rscript1} {params.outdir}
cp {params.rscript2} {params.outdir}
cd {params.outdir}
module load {params.rver}
Rscript limmacall.R '{params.outdir}' '{input.file1}' '{input.file2}' '{params.contrasts}' '{params.refer}' '{params.projectId}' '{params.projDesc}' '{params.gtffile}' '{params.dtype}' '{params.karyobeds}'
"""

rule vennDiagram:
  input:
    join(workpath,"DEG_{dtype}","limma_MDS.png"),
    join(workpath,"DEG_{dtype}","edgeR_prcomp.png"),
    join(workpath,"DEG_{dtype}","DESeq2_PCA.png"),
    limmafile=join(workpath,"DEG_{dtype}","limma_DEG_{group1}-{group2}_all_genes.txt"),
    deseq2file=join(workpath,"DEG_{dtype}","edgeR_DEG_{group1}-{group2}_all_genes.txt"),
    edgeRfile=join(workpath,"DEG_{dtype}","DESeq2_deg_{group1}_vs_{group2}.txt"),
  output: 
    join(workpath,"DEG_{dtype}","limma_edgeR_DESeq2_{group1}-{group2}_vennDiagram.png")
  params:
    rname='pl:vennDiagram',
    batch='--mem=12g --time=2:00:00',
    outdir=join(workpath,"DEG_{dtype}"),
    dtype="{dtype}",
    rver=config['bin'][pfamily]['tool_versions']['RVER'],
    rscript=join(workpath,"Scripts","limma_edgeR_DESeq2_venn.R"),
  shell: """
cp {params.rscript} {params.outdir}
cd {params.outdir}
module load {params.rver}
Rscript {params.rscript} --contrast '{wildcards.group1}-{wildcards.group2}' --limma '{input.limmafile}' --edgeR '{input.edgeRfile}' --DESeq2 '{input.deseq2file}'
"""

rule pca:
  input: 
    file1=join(workpath,star_dir,"sampletable.txt"),
    file2=join(workpath,"DEG_{dtype}","RawCountFile_{dtype}_filtered.txt"),
  output: 
    join(workpath,"DEG_{dtype}","PcaReport.html")
  params: 
    rname='pl:pca',
    batch='--mem=24g --time=10:00:00',
    outdir=join(workpath,"DEG_{dtype}"),
    contrasts=" ".join(config['project']['contrasts']['rcontrasts']),
    dtype="{dtype}",
    projectId=config['project']['id'],
    projDesc=config['project']['description'],
    rver=config['bin'][pfamily]['tool_versions']['RVER'],
    rscript1=join(workpath,"Scripts","pcacall.R"),
    rscript2=join(workpath,"Scripts","PcaReport.Rmd"),
  shell: """
cp {params.rscript1} {params.outdir}
cp {params.rscript2} {params.outdir}
cd {params.outdir}
module load {params.rver}
Rscript pcacall.R '{params.outdir}' '{input.file1}' '{input.file2}' '{params.projectId}' '{params.projDesc}'
"""


# rule salmon:
#   input: bam="{name}.p2.Aligned.toTranscriptome.out.bam"
#   output: "salmonrun/{name}/quant.sf"
#   params: sname="{name}",rname='pl:salmon',batch='--mem=128g --cpus-perptask=8 --time=10:00:00',dir=config['project']['workpath'],rsemref=config['references'][pfamily]['SALMONREF'],libtype={0:'U',1:'SF',2:'SR'}.get(config['bin'][pfamily]['STRANDED'])
#   shell: "mkdir -p {params.dir}/salmonrun; module load salmon/0.6.0; salmon quant -t {params.rsemref} -l I{params.libtype} -a {input.bam} -o {params.dir}/salmonrun/{params.sname} --numBootstraps 30;"

# rule sleuth:
#   input: samtab = "sampletable.txt", bam=expand("salmonrun/{name}/quant.sf", name=samples)
#   output: "salmonrun/sleuth_completed.txt"
#   params: rname='pl:sleuth',batch='--mem=128g --cpus-per-task=8 --time=10:00:00',dir=config['project']['workpath'],pipeRlib=config['bin'][pfamily]['PIPERLIB'],contrasts=" ".join(config['project']['contrasts']['rcontrasts']),species=config['project']['annotation']
#   shell: "module load R/3.5; Rscript Scripts/sleuth.R '{params.dir}' '{params.pipeRlib}' '{input.samtab}' '{params.contrasts}' '{params.species}'"
# 
rule init_deg:
	input:
		samtab=join(workpath,star_dir,"sampletable.txt"),
		rawcountstab=join(workpath,"DEG_ALL","RawCountFile_RSEM_genes.txt"),
	output:
		folders=expand(join(workpath,"DEG_{deg_dir}"),deg_dir=contrasts_w_cpm_cutoff_list),
		samtabs=expand(join(workpath,"DEG_{deg_dir}","sampletable.txt"),deg_dir=contrasts_w_cpm_cutoff_list),
		filteredcountstabs=expand(join(workpath,"DEG_{deg_dir}","RawCountFile_RSEM_genes_filtered.txt"),deg_dir=contrasts_w_cpm_cutoff_list),
	params:
		rname="pl:init_deg_dir",
		rscript=join(workpath,"Scripts","filtersamples.R"),
		rver=config['bin'][pfamily]['tool_versions']['RVER'],
		pyscript=join(workpath,"Scripts","filtersamples.R"),
	shell:"""
for f in {output.folders};do
if [ ! -d $f ]; then mkdir $f ;fi
python {params.pyscript} $f {params.rver} {params.rscript} {input.samtab} {input.rawcountstab} $f/sampletable.txt $f/RawCountFile_RSEM_genes_filtered.txt
done
"""

rule EBSeq_isoform:
  input: 
    samtab=join(workpath,star_dir,"sampletable.txt"),
    samtab2s=expand(join(workpath,"DEG_{deg_dir}","sampletable.txt"),deg_dir=contrasts_w_cpm_cutoff_list),
    igfiles=expand(join(workpath,"DEG_ALL","{name}.RSEM.isoforms.results"),name=samples),
    outdirs=expand(join(workpath,"DEG_{deg_dir}"),deg_dir=contrasts_w_cpm_cutoff_list),
  output: 
    join(workpath,"DEG_{deg_dir}","EBSeq_isoform_completed.txt"),
  params:
  	rname='pl:EBSeq_isoform',
  	batch='--mem=128g --cpus-per-task=8 --time=10:00:00',
  	workpath=workpath,
  	outdir=join(workpath,"DEG_{deg_dir}"),
  	contrasts=" ".join(config['project']['contrasts']['rcontrasts']),
  	rsemref=config['references'][pfamily]['RSEMREF'],
  	rsem=config['bin'][pfamily]['RSEM'],
  	annotate=config['references'][pfamily]['ANNOTATEISOFORMS'],
  	pythonver=config['bin'][pfamily]['tool_versions']['PYTHONVER'],
  	rsemver=config['bin'][pfamily]['tool_versions']['RSEMVER'],
  	script1=join(workpath,"Scripts","EBSeq.py"),
  shell: """
cd {params.outdir}
module load {params.pythonver}
module load {params.rsemver}
python {params.script1} '{params.outdir}' '{input.igfiles}' '{input.samtab}' '{params.contrasts}' '{params.rsemref}' 'isoform' '{params.annotate}'
"""

rule EBSeq_gene:
  input: 
    samtab=join(workpath,star_dir,"sampletable.txt"),
    igfiles=expand(join(workpath,"DEG_ALL","{name}.RSEM.genes.results"), name=samples),
  output: 
    join(workpath,"DEG_{deg_dir}","EBSeq_gene_completed.txt"),
  params:
    rname='pl:EBSeq_gene',
    batch='--mem=128g --cpus-per-task=8 --time=10:00:00',
    outdir=join(workpath,"DEG_{deg_dir}"),
    contrasts=" ".join(config['project']['contrasts']['rcontrasts']),
    rsemref=config['references'][pfamily]['RSEMREF'],
    rsem=config['bin'][pfamily]['RSEM'],
    annotate=config['references'][pfamily]['ANNOTATE'],
    pythonver=config['bin'][pfamily]['tool_versions']['PYTHONVER'],
    rsemver=config['bin'][pfamily]['tool_versions']['RSEMVER'],
    script1=join(workpath,"Scripts","EBSeq.py"),
  shell: """
cd {params.outdir}
module load {params.pythonver}
module load {params.rsemver}
python {params.script1} '{params.outdir}' '{input.igfiles}' '{input.samtab}' '{params.contrasts}' '{params.rsemref}' 'gene' '{params.annotate}'
"""

