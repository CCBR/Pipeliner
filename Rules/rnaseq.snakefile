from snakemake.utils import R
from os.path import join


configfile: "run.json"

samples=config['project']['groups']['rsamps']

from os import listdir
mypath = config['project']['workpath']
fR1 = [f for f in listdir(mypath) if f.find("R1") > 0]
fR2 = [f for f in listdir(mypath) if f.find("R2") > 0]
pe = ""
se = ""
if len(fR1) >= 2 and len(fR2) >= 2:
   pe="yes"
elif len(fR1) >= 2:
   se="yes"
else:
   pass

workpath = config['project']['workpath']
print("Workpath is %s"%(workpath))
if pe=="yes":
  print("Pair-end raw files found.")
elif se=="yes":
  print("Single-end raw files found.")
trim_dir='trim'
star_dir="STAR_files"
bams_dir="bams"
log_dir="logfiles"
rseqc_dir="RSeQC"
rsemi_dir="DEG_RSEM_isoforms"
rsemg_dir="DEG_RSEM_genes"
salmon_dir="Salmon"
subreadg_dir="DEG_Subread_genes"
subreadj_dir="DEG_Subread_junctions"
subreadgj_dir="DEG_Subread_genejunctions"




if config['project']['DEG'] == "yes" and config['project']['TRIM'] == "yes":
  rule all:
     params: 
      batch='--time=168:00:00',
      # input: "STAR_QC",
      # "Reports/multiqc_report.html",
     input: 
      join(workpath,"Reports","multiqc_report.html"),
      join(workpath,star_dir,"sampletable.txt"),
      
      #EBSeq
      
      # join(workpath,rsemi_dir,"EBSeq_isoform_completed.txt"),
      # join(workpath,rsemg_dir,"EBSeq_gene_completed.txt"),
      
      #Salmon/sleuth
      
      # join(workpath,salmon_dir,"sleuth_completed.txt",
      
      #Subread-genes
      
      join(workpath,subreadg_dir,"RawCountFile_Subread_genes_filtered.txt"),
      join(workpath,subreadg_dir,"DESeq2_PCA.png"),
      join(workpath,subreadg_dir,"edgeR_prcomp.png"),
      join(workpath,subreadg_dir,"limma_MDS.png"),
      join(workpath,subreadg_dir,"PcaReport.html"),
      
      #Subread-junctions
      
      join(workpath,subreadj_dir,"RawCountFile_Subread_junctions_filtered.txt"),
      join(workpath,subreadj_dir,"DESeq2_PCA.png"),
      join(workpath,subreadj_dir,"edgeR_prcomp.png"),
      join(workpath,subreadj_dir,"limma_MDS.png"),
      join(workpath,subreadj_dir,"PcaReport.html"),
      
      #Subread-gene-junctions

      join(workpath,subreadgj_dir,"RawCountFile_Subread_genejunctions_filtered.txt"),
      join(workpath,subreadgj_dir,"DESeq2_PCA.png"),
      join(workpath,subreadgj_dir,"edgeR_prcomp.png"),
      join(workpath,subreadgj_dir,"limma_MDS.png"),
      join(workpath,subreadgj_dir,"PcaReport.html"),

      #RSEM
      
      join(workpath,rsemg_dir,"RawCountFile_RSEM_genes_filtered.txt"), 
      join(workpath,rsemg_dir,"DESeq2_PCA.png"),
      join(workpath,rsemg_dir,"edgeR_prcomp.png"),
      join(workpath,rsemg_dir,"limma_MDS.png"),
      join(workpath,rsemg_dir,"PcaReport.html"),
      expand(join(workpath,rsemg_dir,"{name}.RSEM.genes.results"),name=samples),
      join(workpath,rsemg_dir,"RSEM.genes.FPKM.all_samples.txt"),
      join(workpath,rsemi_dir,"RSEM.isoforms.FPKM.all_samples.txt"),

      expand(join(workpath,star_dir,"{name}.star.count.overlap.txt"),name=samples),
      join(workpath,star_dir,"RawCountFileOverlap.txt"),
      join(workpath,star_dir,"RawCountFileStar.txt"),
            
        
elif config['project']['DEG'] == "no" and config['project']['TRIM'] == "yes":
  rule all:
     params: batch='--time=168:00:00'
#     input: "STAR_QC","Reports/multiqc_report.html",
     input: "Reports/multiqc_report.html",
#            expand("{name}.RnaSeqMetrics.txt",name=samples),
#            "postTrimQC",
            "sampletable.txt",
            "RawCountFile_genes_filtered.txt",
            "RawCountFile_junctions_filtered.txt",
            "RawCountFile_genejunctions_filtered.txt",
            expand("{name}.star.count.overlap.txt",name=samples),
            "RawCountFileOverlap.txt","RawCountFileStar.txt",
            expand("{name}.rsem.genes.results",name=samples),
            "DEG_genes/PcaReport.html","DEG_junctions/PcaReport.html","DEG_genejunctions/PcaReport.html","DEG_rsemgenes/PcaReport.html"

elif config['project']['DEG'] == "yes" and config['project']['TRIM'] == "no":
  rule all:
#     input: "STAR_QC","Reports/multiqc_report.html",
     input: "Reports/multiqc_report.html",
            "RSEMDEG_isoform/ebseq_isoform_completed.txt",
            "RSEMDEG_gene/ebseq_gene_completed.txt",
            "salmonrun/sleuth_completed.txt",
            expand("{name}.RnaSeqMetrics.txt",name=samples),
            "sampletable.txt",
            "DEG_genes/deseq2_pca.png",
            "DEG_genes/edgeR_prcomp.png",
            "RawCountFile_genes_filtered.txt",
            "DEG_genes/Limma_MDS.png",
            "DEG_junctions/deseq2_pca.png",
            "DEG_junctions/edgeR_prcomp.png",
            "RawCountFile_junctions_filtered.txt",
            "DEG_junctions/Limma_MDS.png",
            "DEG_genejunctions/deseq2_pca.png",
            "DEG_genejunctions/edgeR_prcomp.png",
            "RawCountFile_genejunctions_filtered.txt",
            "DEG_genejunctions/Limma_MDS.png",
            expand("{name}.star.count.overlap.txt",name=samples),
            "RawCountFileOverlap.txt","RawCountFileStar.txt",
            expand("{name}.rsem.genes.results",name=samples),
            "DEG_genes/PcaReport.htm","DEG_junctions/PcaReport.html","DEG_genejunctions/PcaReport.html"
            
     params: batch='--time=168:00:00'
     #input: "files_to_rnaseqc.txt","STAR_QC","RawCountFile_filtered.txt","sampletable.txt","deseq2_pca.png","edgeR_prcomp.png","Limma_MDS.png"
else:
  rule all:
     params: batch='--time=168:00:00'
#     input: "STAR_QC","Reports/multiqc_report.html",
     input: "Reports/multiqc_report.html",
            expand("{name}.RnaSeqMetrics.txt",name=samples),
            "RawCountFile_genes_filtered.txt",
            "RawCountFile_genejunctions_filtered.txt",
            expand("{name}.star.count.overlap.txt",name=samples),
            "RawCountFileOverlap.txt","RawCountFileStar.txt",
            expand("{name}.rsem.genes.results",name=samples),
            "DEG_genes/PcaReport.htm","DEG_junctions/PcaReport.html","DEG_genejunctions/PcaReport.html"

########################################################

##############################################

rule get_strandness:
  input: 
    groupsfile=join(workpath,"groups.tab"),
  output: 
    outfile=join(workpath,log_dir,"strandness.txt"),
    outdir=join(workpath,log_dir)
  params: 
    rname='pl:get_strandness',
    pythonver=config['bin'][pfamily]['tool_versions']['PYTHONVER'],
    pythonscript=join(workpath,"Scripts","get_strandness.py")
  run:
    import os
    os.chdir(output.outdir)
    os.system("module load "+params.pythonver+";python "+params.pythonscript+" "+input.groupsfile+" > "+output.outfile)
    strandfile=open(output.outfile,'r')
    strandness=strandfile.readline().strip()
    strandfile.close()
    A=open(join(workpath,"run.json"),'r')
    a=eval(A.read())
    A.close()
    config=dict(a.items())
    config['project']['STRANDED']=strandness
    with open(join(workpath,'run.json'),'w') as F:
      json.dump(config, F, sort_keys = True, indent = 4,ensure_ascii=False)
    F.close()
   
rule samplecondition:
   input: 
    files=expand(join(workpath,star_dir,"{name}.star.count.txt"), name=samples)
   output: 
    out1=join(workpath,star_dir,"sampletable.txt")
   params: 
    rname='pl:samplecondition',
    batch='--mem=4g --time=10:00:00', 
    groups=config['project']['groups']['rgroups'],
    labels=config['project']['groups']['rlabels'],
    gtffile=config['references'][pfamily]['GTFFILE']
   run:
        with open(output.out1, "w") as out:
            out.write("sampleName\tfileName\tcondition\tlabel\n")
            i=0
            for f in input.files:
                out.write("%s\t"  % f)
                out.write("%s\t"  % f)
                out.write("%s\t" % params.groups[i])
                out.write("%s\n" % params.labels[i])                
                i=i+1
            out.close()

if pe=="yes":

   rule rsem:
      input: 
        file1=join(workpath,bams_dir,"{name}.p2.Aligned.toTranscriptome.out.bam"),
      output: 
        out1=join(workpath,rsemg_dir,"{name}.RSEM.genes.results"),
        out2=join(workpath,rsemg_dir,"{name}.RSEM.isoforms.results"),
      params:
        rname='pl:rsem',
        prefix="{name}.RSEM",
        outdir=join(workpath,rsemg_dir),
        outdir2=join(workpath,rsemi_dir),
        batch='--cpus-per-task=16 --mem=32g --time=24:00:00',
        rsemref=config['references'][pfamily]['RSEMREF'],
        rsemver=config['bin'][pfamily]['tool_versions']['RSEMVER'],
        pythonver=config['bin'][pfamily]['tool_versions']['PYTHONVER'],
        annotate=config['references'][pfamily]['ANNOTATE'],
        pythonscript=join(workpath,"Scripts","merge_rsem_results.py"),
      threads: 16
      shell: """
cd {params.outdir}
module load {params.rsemver}
rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  --bam --paired-end -p {threads}  {input.file1} {params.rsemref} {params.prefix} --time --temporary-folder /lscratch/$SLURM_JOBID --keep-intermediate-files
mv {params.outdir}/{params.prefix}.isoforms.results {params.outdir2}
"""

if se=="yes":

   rule rsem:
      input:
        file1=join(workpath,bams_dir,"{name}.p2.Aligned.toTranscriptome.out.bam"),
      output:
        out1=join(workpath,rsemg_dir,"{name}.RSEM.genes.results"),
        out2=join(workpath,rsemi_dir,"{name}.RSEM.isoforms.results"),
      params: 
        rname='pl:rsem',
        prefix="{name}.RSEM",
        outdir=join(workpath,rsemg_dir),
        outdir2=join(workpath,rsemi_dir),
        batch='--cpus-per-task=16 --mem=32g --time=24:00:00',
        rsemref=config['references'][pfamily]['RSEMREF'],
        rsemver=config['bin'][pfamily]['tool_versions']['RSEMVER'],
        pythonver=config['bin'][pfamily]['tool_versions']['PYTHONVER'],
        annotate=config['references'][pfamily]['ANNOTATE'],
        pythonscript=join(workpath,"Scripts","merge_rsem_results.py"),
      threads: 16
      shell: """
cd {params.outdir}
module load {params.rsemver}
rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  --bam -p {threads}  {input.file1} {params.rsemref} {params.prefix} --time --temporary-folder /lscratch/$SLURM_JOBID --keep-intermediate-files
mv {params.outdir}/{params.prefix}.isoforms.results {params.outdir2}
"""

rule rsem_merge:
  input:
    join(workpath,rsemg_dir),
    join(workpath,rsemi_dir),
  output:
    fpkm1=join(workpath,rsemg_dir,"RSEM.genes.FPKM.all_samples.txt"),
    fpkm2=join(workpath,rsemi_dir,"RSEM.isoforms.FPKM.all_samples.txt"),
  params: 
    rname='pl:rsem_merge',
    pythonver=config['bin'][pfamily]['tool_versions']['PYTHONVER'],
    annotate=config['references'][pfamily]['ANNOTATE'],
    pythonscript=join(workpath,"Scripts","merge_rsem_results.py"),
  threads: 16
  shell: """
module load {params.pythonver}
python {params.pythonscript} {params.annotate} {rsemg_dir} {rsemi_dir}
"""




rule rsemcounts:
   input:
    files=expand(join(workpath,rsemg_dir,"{name}.RSEM.genes.results"), name=samples),
   output: 
    join(workpath,rsemg_dir,"RawCountFile_RSEM_genes_filtered.txt"),
   params: 
    rname='pl:rsemcounts',
    batch='--mem=8g --time=10:00:00',
    outdir=join(workpath,rsemg_dir),
    mincount=config['project']['MINCOUNTGENES'],
    minsamples=config['project']['MINSAMPLES'],
    annotate=config['references'][pfamily]['ANNOTATE'],
    rver=config['bin'][pfamily]['tool_versions']['RVER'],
    rscript=join(workpath,"Scripts","rsemcounts.R")
   shell: """
cd {params.outdir}
module load {params.rver}
Rscript {params.rscript} '{params.outdir}' '{input.files}' '{params.mincount}' '{params.minsamples}' '{params.annotate}'
"""


rule subread:
   input:
    file1=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
    file2=join(workpath,log_dir,"strandness.txt")
   output:
    out=join(workpath,star_dir,"{name}.star.count.info.txt"),
    res=join(workpath,star_dir,"{name}.star.count.txt"),
   params:
    rname='pl:subread',
    batch='--time=4:00:00 --gres=lscratch:800',
    subreadver=config['bin'][pfamily]['tool_versions']['SUBREADVER'],
    gtffile=config['references'][pfamily]['GTFFILE'],
   threads: 16
   shell: """
module load {params.subreadver}
featureCounts -T {threads} -s `cat {input.file2}` -p -t exon -R -g gene_id -a {params.gtffile} --tmpDir /lscratch/$SLURM_JOBID  -o {output.out}  {input.file1}
sed '1d' {output.out} | cut -f1,7 > {output.res}
"""

rule subreadoverlap:
   input: 
    file1=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
    file2=join(workpath,log_dir,"strandness.txt")
   output:
    out=join(workpath,star_dir,"{name}.star.count.info.overlap.txt"),
    res=join(workpath,star_dir,"{name}.star.count.overlap.txt")
   params: 
    rname='pl:subreadoverlap',
    batch='--cpus-per-task=16 --mem=24g --time=48:00:00 --gres=lscratch:800',
    subreadver=config['bin'][pfamily]['tool_versions']['SUBREADVER'],
    gtffile=config['references'][pfamily]['GTFFILE']
   threads: 16
   shell: """
module load {params.subreadver}
featureCounts -T {threads} -s `cat {input.file2}` -p -t exon -R -O -g gene_id -a {params.gtffile} --tmpDir /lscratch/$SLURM_JOBID  -o {output.out}  {input.file1}
sed '1d' {output.out} | cut -f1,7 > {output.res}
"""

rule genecounts: 
   input:
    file1=expand(join(workpath,star_dir,"{name}.star.count.txt"), name=samples),
    file2=join(workpath,star_dir,"sampletable.txt")
   output:
    join(workpath,subreadg_dir,"RawCountFile_Subread_genes_filtered.txt")
   params: 
    rname='pl:genecounts',
    batch='--mem=8g --time=10:00:00',
    outdir=join(workpath,subreadg_dir),
    mincount=config['project']['MINCOUNTGENES'],
    minsamples=config['project']['MINSAMPLES'],
    annotate=config['references'][pfamily]['ANNOTATE'],
    rver=config['bin'][pfamily]['tool_versions']['RVER'],
    rscript=join(workpath,"Scripts","genecounts.R"),
   shell: """
module load {params.rver}
Rscript {params.rscript} '{params.outdir}' '{input.file1}' '{params.mincount}' '{params.minsamples}' '{params.annotate}' '{input.file2}'
"""

rule junctioncounts: 
   input: 
    files=expand(join(workpath,star_dir,"{name}.p2.SJ.out.tab"), name=samples)
   output:
    join(workpath,subreadj_dir,"RawCountFile_Subread_junctions_filtered.txt")
   params:
    rname='pl:junctioncounts',
    batch='--mem=8g --time=10:00:00',
    outdir=join(workpath,subreadj_dir),
    mincount=config['project']['MINCOUNTJUNCTIONS'],
    minsamples=config['project']['MINSAMPLES'],
    rver=config['bin'][pfamily]['tool_versions']['RVER'],
    rscript=join(workpath,"Scripts","junctioncounts.R"),
   shell: """
module load {params.rver}
Rscript {params.rscript} '{params.outdir}' '{input.files}' '{params.mincount}' '{params.minsamples}'
"""

rule genejunctioncounts: 
   input:
    files=expand(join(workpath,star_dir,"{name}.p2.SJ.out.tab"), name=samples)
   output: 
    join(workpath,subreadgj_dir,"RawCountFile_Subread_genejunctions_filtered.txt")
   params:
    rname='pl:genejunctions',
    batch='--mem=8g --time=10:00:00',
    outdir=join(workpath,subreadgj_dir),
    geneinfo=config['references'][pfamily]['GENEINFO'],
    mincount=config['project']['MINCOUNTGENEJUNCTIONS'],
    minsamples=config['project']['MINSAMPLES'],
    rver=config['bin'][pfamily]['tool_versions']['RVER'],
    rscript=join(workpath,"Scripts","genejunctioncounts.R")
   shell: """
module load {params.rver}
Rscript {params.rscript} '{params.outdir}' '{input.files}' '{params.geneinfo}' '{params.mincount}' '{params.minsamples}'
"""

rule joincounts:
   input: 
    files=expand(join(workpath,star_dir,"{name}.star.count.overlap.txt"), name=samples),
    files2=expand(join(workpath,star_dir,"{name}.p2.ReadsPerGene.out.tab"), name=samples)
   output: 
    out1=join(workpath,star_dir,"RawCountFileOverlap.txt"),
    out2=join(workpath,star_dir,"RawCountFileStar.txt")
   params: 
    rname='pl:junctioncounts',
    batch='--mem=8g --time=10:00:00',
    outdir=join(workpath,star_dir),
    starstrandcol=config['bin'][pfamily]['STARSTRANDCOL'],
    rver=config['bin'][pfamily]['tool_versions']['RVER'],
    rscript=join(workpath,"Scripts","joincounts.R")
   shell: """
module load {params.rver}
Rscript {params.rscript} '{params.outdir}' '{input.files}' '{input.files2}' '{params.starstrandcol}'
"""


rule rnaseq_multiqc:
#    input: "STAR_QC/index.html","STAR_QC/report.html"
    input:
      join(workpath,star_dir,"sampletable.txt")
    output:
      join(workpath,"Reports","multiqc_report.html")
    params:
      rname="pl:multiqc",
      outdir=join(workpath,"Reports"),
      multiqcver=config['bin'][pfamily]['tool_versions']['MULTIQCVER'],
      qcconfig=config['bin'][pfamily]['CONFMULTIQC']
    threads: 1
    shell:  """
module load {params.multiqcver}
cd {params.outdir}
multiqc -f -c {params.qcconfig}  ../
"""

rule deseq2:
  input:
    file1=join(workpath,star_dir,"sampletable.txt"),
    file2=join(workpath,"DEG_{dtype}","RawCountFile_{dtype}_filtered.txt"),
  output:
    join(workpath,"DEG_{dtype}","DESeq2_PCA.png"),
  params:
    rname='pl:deseq2',
    batch='--mem=24g --time=10:00:00',
    outdir=join(workpath,"DEG_{dtype}"),
    contrasts=" ".join(config['project']['contrasts']['rcontrasts']),
    dtype="{dtype}",
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
    file1=join(workpath,star_dir,"sampletable.txt"),
    file2=join(workpath,"DEG_{dtype}","RawCountFile_{dtype}_filtered.txt"),
  output:
    join(workpath,"DEG_{dtype}","edgeR_prcomp.png"),
  params: 
    rname='pl:edgeR',
    batch='--mem=24g --time=10:00:00',
    outdir=join(workpath,"DEG_{dtype}"),
    contrasts=" ".join(config['project']['contrasts']['rcontrasts']),
    dtype="{dtype}",
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
    file1=join(workpath,star_dir,"sampletable.txt"),
    file2=join(workpath,"DEG_{dtype}","RawCountFile_{dtype}_filtered.txt"),
  output: 
    join(workpath,"DEG_{dtype}","limma_MDS.png")
  params:
    rname='pl:limmavoom',
    batch='--mem=24g --time=10:00:00',
    outdir=join(workpath,"DEG_{dtype}"),
    contrasts=" ".join(config['project']['contrasts']['rcontrasts']),
    dtype="{dtype}",
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
#   shell: "module load R/3.4.0_gcc-6.2.0; Rscript Scripts/sleuth.R '{params.dir}' '{params.pipeRlib}' '{input.samtab}' '{params.contrasts}' '{params.species}'"

rule EBSeq_isoform:
  input: 
    samtab=join(workpath,star_dir,"sampletable.txt"),
    igfiles=expand(join(workpath,rsemi_dir,"{name}.RSEM.isoforms.results"), name=samples),
  output: 
    join(workpath,rsemi_dir,"EBSeq_isoform_completed.txt")
  params: 
    rname='pl:EBSeq',
    batch='--mem=128g --cpus-per-task=8 --time=10:00:00',
    outdir=join(workpath,rsemi_dir),
    contrasts=" ".join(config['project']['contrasts']['rcontrasts']),
    rsemref=config['references'][pfamily]['RSEMREF'],
    rsem=config['bin'][pfamily]['RSEM'],
    annotate=config['references'][pfamily]['ANNOTATE'],
    rver=config['bin'][pfamily]['tool_versions']['RVER'],
    rsemver=config['bin'][pfamily]['tool_versions']['RSEMVER'],
    rscript1=join(workpath,"Scripts","ebseq.R"),
  shell: """
cp {params.rscript1} {params.outdir}
cd {params.outdir}
module load {params.rver}
module load {params.rsemver}
Rscript ebseq.R '{params.outdir}' '{input.igfiles}' '{input.samtab}' '{params.contrasts}' '{params.rsemref}' 'isoform' '{params.annotate}'
"""

rule EBSeq_gene:
  input: 
    samtab=join(workpath,star_dir,"sampletable.txt"),
    igfiles=expand(join(workpath,rsemg_dir,"{name}.RSEM.genes.results"), name=samples),
  output: 
    join(workpath,rsemg_dir,"EBSeq_gene_completed.txt")
  params:
    rname='pl:EBSeq',
    batch='--mem=128g --cpus-per-task=8 --time=10:00:00',
    outdir=join(workpath,rsemg_dir),
    contrasts=" ".join(config['project']['contrasts']['rcontrasts']),
    rsemref=config['references'][pfamily]['RSEMREF'],
    rsem=config['bin'][pfamily]['RSEM'],
    annotate=config['references'][pfamily]['ANNOTATE'],
    rver=config['bin'][pfamily]['tool_versions']['RVER'],
    rsemver=config['bin'][pfamily]['tool_versions']['RSEMVER'],
    rscript1=join(workpath,"Scripts","ebseq.R"),
  shell: """
cp {params.rscript1} {params.outdir}
cd {params.outdir}
module load {params.rver}
module load {params.rsemver}
Rscript ebseq.R '{params.outdir}' '{input.igfiles}' '{input.samtab}' '{params.contrasts}' '{params.rsemref}' 'gene' '{params.annotate}'
"""
