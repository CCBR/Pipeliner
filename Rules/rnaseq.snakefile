from snakemake.utils import R

configfile: "run.json"

#samples=config['project']['units']
#samples=config['project']['contrasts']['rsamps']
samples=config['project']['groups']['rsamps']

if config['project']['DEG'] == "yes" and config['project']['TRIM'] == "yes":
  rule all:
     params: batch='--time=168:00:00'
#     input: "STAR_QC","Reports/multiqc_report.html",
     input: "Reports/multiqc_report.html",
            "ebseq_completed.txt",
            "salmonrun/sleuth_completed.txt",
         #   expand("{name}.RnaSeqMetrics.txt",name=samples),
         #   "postTrimQC","sampletable.txt",
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
            "RawCountFileOverlap.txt",
            "RawCountFileStar.txt",expand("{name}.rsem.genes.results",name=samples),
            "DEG_genes/PcaReport.html","DEG_junctions/PcaReport.html","DEG_genejunctions/PcaReport.html"
            
        
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
#            expand("{name}.rsem.genes.results",name=samples),
            "DEG_genes/PcaReport.html","DEG_junctions/PcaReport.html","DEG_genejunctions/PcaReport.html"

elif config['project']['DEG'] == "yes" and config['project']['TRIM'] == "no":
  rule all:
#     input: "STAR_QC","Reports/multiqc_report.html",
     input: "Reports/multiqc_report.html",
            "ebseq_completed.txt",
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
  input: "groups.tab"
  output: "strandness.txt"
  params: rname='pl:get_strandness'
#  shell: "python Scripts/get_strandness.py > strandness.txt"
  run:
    import os
    os.system("python Scripts/get_strandness.py > strandness.txt")
    strandfile=open("strandness.txt",'r')
    strandness=strandfile.readline().strip()
    strandfile.close()
    A=open("run.json",'r')
    a=eval(A.read())
    A.close()
    config=dict(a.items())
    config['project']['STRANDED']=strandness
    with open('run.json','w') as F:
      json.dump(config, F, sort_keys = True, indent = 4,ensure_ascii=False)
    F.close()

rule rsem:
  input: file1= "{name}.p2.Aligned.toTranscriptome.out.bam"
  output: out1="{name}.rsem.genes.results",out2="{name}.rsem.isoforms.results"
  params: rname='pl:rsem',prefix="{name}.rsem",batch='--cpus-per-task=16 --mem=32g --time=24:00:00',rsemref=config['references'][pfamily]['RSEMREF'],rsem=config['bin'][pfamily]['RSEM']
  shell: "{params.rsem}/rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  --bam --paired-end -p 16  {input.file1} {params.rsemref} {params.prefix}"


rule subread:
   input:  file1="{name}.star_rg_added.sorted.dmark.bam", file2="strandness.txt"
   output: out="{name}.star.count.info.txt", res="{name}.star.count.txt"
   params: rname='pl:subread',batch='--time=4:00:00',subreadver=config['bin'][pfamily]['SUBREADVER'],gtffile=config['references'][pfamily]['GTFFILE']
   shell: "module load {params.subreadver}; featureCounts -T 16 -s `cat strandness.txt` -p -t exon -R -g gene_id -a {params.gtffile} -o {output.out}  {input.file1}; sed '1d' {output.out} | cut -f1,7 > {output.res}"

rule subreadoverlap:
   input:  file1="{name}.star_rg_added.sorted.dmark.bam", file2="strandness.txt"
   output: out="{name}.star.count.info.overlap.txt", res="{name}.star.count.overlap.txt"
   params: rname='pl:subreadoverlap',batch='--cpus-per-task=16 --mem=24g --time=48:00:00',subreadver=config['bin'][pfamily]['SUBREADVER'],gtffile=config['references'][pfamily]['GTFFILE']
   shell: "module load {params.subreadver}; featureCounts -T 16 -s `cat strandness.txt` -p -t exon -R -O -g gene_id -a {params.gtffile} -o {output.out}  {input.file1}; sed '1d' {output.out} | cut -f1,7 > {output.res}"


rule genecounts: 
   input: files=expand("{name}.star.count.txt", name=samples)
   output: "RawCountFile_genes_filtered.txt"
   params: rname='pl:genecounts',batch='--mem=8g --time=10:00:00',dir=config['project']['workpath'],mincount=config['project']['MINCOUNTGENES'],minsamples=config['project']['MINSAMPLES'],annotate=config['references'][pfamily]['ANNOTATE']
   shell: "module load R; Rscript Scripts/genecounts.R '{params.dir}' '{input.files}' '{params.mincount}' '{params.minsamples}' '{params.annotate}'"

rule junctioncounts: 
   input: files=expand("{name}.p2.SJ.out.tab", name=samples)
   output: "RawCountFile_junctions_filtered.txt"
   params: rname='pl:junctioncounts',batch='--mem=8g --time=10:00:00',dir=config['project']['workpath'],mincount=config['project']['MINCOUNTJUNCTIONS'],minsamples=config['project']['MINSAMPLES']
   shell: "module load R; Rscript Scripts/junctioncounts.R '{params.dir}' '{input.files}' '{params.mincount}' '{params.minsamples}'"

rule genejunctioncounts: 
   input: files=expand("{name}.p2.SJ.out.tab", name=samples)
   output: "RawCountFile_genejunctions_filtered.txt"
   params: rname='pl:genejunctions',batch='--mem=8g --time=10:00:00',dir=config['project']['workpath'],gtffile=config['references'][pfamily]['GTFFILE'],mincount=config['project']['MINCOUNTGENEJUNCTIONS'],minsamples=config['project']['MINSAMPLES']
   shell: "module load R; Rscript Scripts/genejunctioncounts.R '{params.dir}' '{input.files}' '{params.gtffile}' '{params.mincount}' '{params.minsamples}'"

rule joincounts:
   input: files=expand("{name}.star.count.overlap.txt", name=samples),files2=expand("{name}.p2.ReadsPerGene.out.tab", name=samples)
   output: "RawCountFileOverlap.txt","RawCountFileStar.txt"
   params: rname='pl:junctioncounts',batch='--mem=8g --time=10:00:00',dir=config['project']['workpath'],starstrandcol=config['bin'][pfamily]['STARSTRANDCOL']
   shell: "module load R; Rscript Scripts/joincounts.R '{params.dir}' '{input.files}' '{input.files2}' '{params.starstrandcol}'"


rule rnaseq_multiqc:
#    input: "STAR_QC/index.html","STAR_QC/report.html"
    input: "sampletable.txt"
    output: "Reports/multiqc_report.html"
#    params: rname="pl:multiqc",pythonpath=config['bin'][pfamily]['PYTHONPATH'],multiqc=config['bin'][pfamily]['MULTIQC']
    params: rname="pl:multiqc",multiqc=config['bin'][pfamily]['MULTIQC'],qcconfig=config['bin'][pfamily]['CONFMULTIQC']
    threads: 1
    shell:  """
#            module load multiqc
#               cd Reports && multiqc -f -e featureCounts -e picard ../
#            cd Reports && multiqc -f -e featureCounts  ../ 
            module load {params.multiqc}
            cd Reports && multiqc -f -c {params.qcconfig}  ../
            """

   
rule samplecondition:
   input: files=expand("{name}.star.count.txt", name=samples)
   output: out1= "sampletable.txt"
   params: rname='pl:samplecondition',batch='--mem=4g --time=10:00:00', groups=config['project']['groups']['rgroups'], labels=config['project']['groups']['rlabels']
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

rule deseq2:
  input: file1="sampletable.txt", file2="RawCountFile{dtype}_filtered.txt"
  ## input: "sampletable.txt"
  output: "DEG{dtype}/deseq2_pca.png"
  params: rname='pl:deseq2',batch='--mem=24g --time=10:00:00',dir=config['project']['workpath'],annotate=config['references'][pfamily]['ANNOTATE'],contrasts=" ".join(config['project']['contrasts']['rcontrasts']), dtype="{dtype}"
##  shell: "mkdir -p DEG{params.dtype}; module load R; Rscript Scripts/deseq2.R '{params.dir}/DEG{params.dtype}/' '../{input.file1}' '../{input.file2}' '{params.annotate}' '{params.contrasts}'"
  shell: "mkdir -p DEG{params.dtype}; cp Scripts/Deseq2Report.Rmd {params.dir}/DEG{params.dtype}/; module load R; Rscript Scripts/deseq2call.R '{params.dir}/DEG{params.dtype}/' '../{input.file1}' '../{input.file2}' '{params.contrasts}'"

rule edgeR:
  input: file1="sampletable.txt", file2="RawCountFile{dtype}_filtered.txt"
  output: "DEG{dtype}/edgeR_prcomp.png"
  params: rname='pl:edgeR',batch='--mem=24g --time=10:00:00', dir=config['project']['workpath'],annotate=config['references'][pfamily]['ANNOTATE'],contrasts=" ".join(config['project']['contrasts']['rcontrasts']), dtype="{dtype}"
##  shell: "mkdir -p DEG{params.dtype}; module load R; Rscript Scripts/edgeR.R '{params.dir}/DEG{params.dtype}/' '../{input.file1}' '../{input.file2}' '{params.annotate}' '{params.contrasts}'"
  shell: "mkdir -p DEG{params.dtype}; cp Scripts/EdgerReport.Rmd {params.dir}/DEG{params.dtype}/; module load R; Rscript Scripts/edgeRcall.R '{params.dir}/DEG{params.dtype}/' '../{input.file1}' '../{input.file2}' '{params.contrasts}'"

rule limmavoom:
  input: file1="sampletable.txt", file2="RawCountFile{dtype}_filtered.txt"
  output: "DEG{dtype}/Limma_MDS.png"
  params: rname='pl:limmavoom',batch='--mem=24g --time=10:00:00',dir=config['project']['workpath'],contrasts=" ".join(config['project']['contrasts']['rcontrasts']), dtype="{dtype}"
##  shell: "mkdir -p DEG{params.dtype}; module load R; Rscript Scripts/limmavoom.R '{params.dir}/DEG{params.dtype}/' '../{input.file1}' '../{input.file2}' '{params.contrasts}'"
  shell: "mkdir -p DEG{params.dtype}; cp Scripts/LimmaReport.Rmd {params.dir}/DEG{params.dtype}/; module load R; Rscript Scripts/limmacall.R '{params.dir}/DEG{params.dtype}/' '../{input.file1}' '../{input.file2}' '{params.contrasts}'"

rule pca:
  input: file1="sampletable.txt", file2="RawCountFile{dtype}_filtered.txt"
  output: "DEG{dtype}/PcaReport.html"
  params: rname='pl:pca',batch='--mem=24g --time=10:00:00',dir=config['project']['workpath'],contrasts=" ".join(config['project']['contrasts']['rcontrasts']), dtype="{dtype}"
  shell: "mkdir -p DEG{params.dtype}; cp Scripts/PcaReport.Rmd {params.dir}/DEG{params.dtype}/; module load R; Rscript Scripts/pcacall.R '{params.dir}/DEG{params.dtype}/' '../{input.file1}' '../{input.file2}' '{params.contrasts}'"


rule salmon:
  input: bam="{name}.p2.Aligned.toTranscriptome.out.bam"
  output: "salmonrun/{name}/quant.sf"
  params: sname="{name}",rname='pl:salmon',batch='--mem=128g --cpus-perptask=8 --time=10:00:00',dir=config['project']['workpath'],rsemref=config['references'][pfamily]['SALMONREF'],libtype={0:'U',1:'SF',2:'SR'}.get(config['bin'][pfamily]['STRANDED'])
  shell: "mkdir -p {params.dir}/salmonrun; module load salmon/0.6.0; salmon quant -t {params.rsemref} -l I{params.libtype} -a {input.bam} -o {params.dir}/salmonrun/{params.sname} --numBootstraps 30;"

rule sleuth:
  input: samtab = "sampletable.txt", bam=expand("salmonrun/{name}/quant.sf", name=samples)
  output: "salmonrun/sleuth_completed.txt"
  params: rname='pl:sleuth',batch='--mem=128g --cpus-per-task=8 --time=10:00:00',dir=config['project']['workpath'],pipeRlib=config['bin'][pfamily]['PIPERLIB'],contrasts=" ".join(config['project']['contrasts']['rcontrasts']),species=config['project']['annotation']
  shell: "module load R; Rscript Scripts/sleuth.R '{params.dir}' '{params.pipeRlib}' '{input.samtab}' '{params.contrasts}' '{params.species}'"

rule EBSeq:
  input: samtab = "sampletable.txt", isoforms=expand("{name}.rsem.isoforms.results", name=samples)
  output: "ebseq_completed.txt"
  params: rname='pl:EBSeq',batch='--mem=128g --cpus-per-task=8 --time=10:00:00',dir=config['project']['workpath'],contrasts=" ".join(config['project']['contrasts']['rcontrasts']),rsemref=config['references'][pfamily]['RSEMREF'],rsem=config['bin'][pfamily]['RSEM']
  shell: "module load R; Rscript Scripts/ebseq.R '{params.dir}' '{input.isoforms}' '{input.samtab}' '{params.contrasts}' '{params.rsemref}' '{params.rsem}'"
