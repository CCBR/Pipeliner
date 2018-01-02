from snakemake.utils import R

configfile: "run.json"

rule trimmomatic_pe:
   input: file1= config['project']['workpath']+"/{name}.R1."+config['project']['filetype'],file2=config['project']['workpath']+"/{name}.R2."+config['project']['filetype']
   output: out11=temp("trim/{name}_R1_001_trim_paired.fastq.gz"),out12=temp("trim/{name}_R1_001_trim_unpaired.fastq.gz"),out21=temp("trim/{name}_R2_001_trim_paired.fastq.gz"),out22=temp("trim/{name}_R2_001_trim_unpaired.fastq.gz"),err="QC/{name}_run_trimmomatic.err"
   params: rname='pl:trimmomatic_pe',batch='--cpus-per-task=32 --mem=110g --time=48:00:00',trimmomaticver=config['bin'][pfamily]['TRIMMOMATICVER'],fastawithadaptersetc=config['references'][pfamily]['FASTAWITHADAPTERSETC'],seedmismatches=config['bin'][pfamily]['SEEDMISMATCHES'],palindromeclipthreshold=config['bin'][pfamily]['PALINDROMECLIPTHRESHOLD'],simpleclipthreshold=config['bin'][pfamily]['SIMPLECLIPTHRESHOLD'],leadingquality=config['bin'][pfamily]['LEADINGQUALITY'],trailingquality=config['bin'][pfamily]['TRAILINGQUALITY'],windowsize=config['bin'][pfamily]['WINDOWSIZE'],windowquality=config['bin'][pfamily]['WINDOWQUALITY'],targetlength=config['bin'][pfamily]['TARGETLENGTH'],strictness=config['bin'][pfamily]['STRICTNESS'],minlen=config['bin'][pfamily]['MINLEN'],headcroplength=config['bin'][pfamily]['HEADCROPLENGTH']
   threads:32
   shell:"module load {params.trimmomaticver}; java -classpath $TRIMMOJAR   org.usadellab.trimmomatic.TrimmomaticPE -threads {threads} {input.file1} {input.file2} {output.out11} {output.out12} {output.out21} {output.out22} ILLUMINACLIP:{params.fastawithadaptersetc}:{params.seedmismatches}:{params.palindromeclipthreshold}:{params.simpleclipthreshold} MINLEN:{params.minlen} 2> {output.err}"

rule star1p:
   input: file1= "trim/{name}_R1_001_trim_paired.fastq.gz",file2="trim/{name}_R2_001_trim_paired.fastq.gz",qcdir="QC",length="QC/{name}_readlength.txt"
   output: out1= temp("{name}.SJ.out.tab"),out3= temp("{name}.Aligned.out.bam")
   params: rname='pl:star1p',prefix="{name}",batch='--cpus-per-task=32 --mem=110g --time=48:00:00',starver=config['bin'][pfamily]['STARVER'],stardir=config['references']['rnaseq']['STARDIR'],filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],filtertype=config['bin'][pfamily]['FILTERTYPE'],filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],adapter1=config['bin'][pfamily]['ADAPTER1'],adapter2=config['bin'][pfamily]['ADAPTER2'],
   threads: 32
   run:
     import glob,json
     rl=int(open(input.length).readlines()[0].strip())-1
     dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.stardir+'*/',recursive=False))))
     bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
     cmd="module load "+params.starver+"; STAR --genomeDir "+params.stardir+str(bestdbrl)+" --outFilterIntronMotifs "+params.filterintronmotifs+" --outSAMstrandField "+params.samstrandfield+"  --outFilterType "+params.filtertype+" --outFilterMultimapNmax "+str(params.filtermultimapnmax)+" --alignSJoverhangMin "+str(params.alignsjoverhangmin)+" --alignSJDBoverhangMin "+str(params.alignsjdboverhangmin)+"  --outFilterMismatchNmax "+str(params.filtermismatchnmax)+" --outFilterMismatchNoverLmax "+str(params.filtermismatchnoverlmax)+"  --alignIntronMin "+str(params.alignintronmin)+" --alignIntronMax "+str(params.alignintronmax)+" --alignMatesGapMax "+str(params.alignmatesgapmax)+" clip3pAdapterSeq "+params.adapter1+" "+params.adapter2+" --readFilesIn "+input.file1+" "+input.file2+" --readFilesCommand zcat --runThreadN "+str(threads)+" --outFileNamePrefix "+params.prefix+". --outSAMtype BAM Unsorted"
      # --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40"
     shell(cmd)
     # add sjdboverhand and stardir back to run.json
     A=open("run.json",'r')
     a=eval(A.read())
     A.close()
     config=dict(a.items())
     config['project']['SJDBOVERHANG']=str(bestdbrl)
     config['project']["STARDIR"]= config['references'][pfamily]['STARDIR']+str(bestdbrl)
     with open('run.json','w') as F:
       json.dump(config, F, sort_keys = True, indent = 4,ensure_ascii=False)
     F.close()

rule star2p:
   input: file1= "trim/{name}_R1_001_trim_paired.fastq.gz",file2="trim/{name}_R2_001_trim_paired.fastq.gz",tab=expand("{name}.SJ.out.tab",name=samples),qcdir="QC",length="QC/{name}_readlength.txt"
   output: out1=temp("{name}.p2.Aligned.sortedByCoord.out.bam"),out2=temp("{name}.p2.ReadsPerGene.out.tab"),out3=temp("{name}.p2.Aligned.toTranscriptome.out.bam"),out4=temp("{name}.p2.SJ.out.tab"),out5=temp("{name}.p2.Log.final.out")
   params: rname='pl:star2p',prefix="{name}.p2",batch='--cpus-per-task=32 --mem=110g --time=48:00:00',starver=config['bin'][pfamily]['STARVER'],filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],filtertype=config['bin'][pfamily]['FILTERTYPE'],filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],adapter1=config['bin'][pfamily]['ADAPTER1'],adapter2=config['bin'][pfamily]['ADAPTER2'],outsamunmapped=config['bin'][pfamily]['OUTSAMUNMAPPED'],wigtype=config['bin'][pfamily]['WIGTYPE'],wigstrand=config['bin'][pfamily]['WIGSTRAND'], gtffile=config['references'][pfamily]['GTFFILE'], nbjuncs=config['bin'][pfamily]['NBJUNCS'],stardir=config['references']['rnaseq']['STARDIR']
   threads:32
   run:
     import glob
     rl=int(open(input.length).readlines()[0].strip())-1
     dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.stardir+'*/',recursive=False))))
     bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
     cmd="module load "+params.starver+"; STAR --genomeDir "+params.stardir+str(bestdbrl)+" --outFilterIntronMotifs "+params.filterintronmotifs+" --outSAMstrandField "+params.samstrandfield+"  --outFilterType "+params.filtertype+" --outFilterMultimapNmax "+str(params.filtermultimapnmax)+" --alignSJoverhangMin "+str(params.alignsjoverhangmin)+" --alignSJDBoverhangMin "+str(params.alignsjdboverhangmin)+"  --outFilterMismatchNmax "+str(params.filtermismatchnmax)+" --outFilterMismatchNoverLmax "+str(params.filtermismatchnoverlmax)+"  --alignIntronMin "+str(params.alignintronmin)+" --alignIntronMax "+str(params.alignintronmax)+" --alignMatesGapMax "+str(params.alignmatesgapmax)+"  --clip3pAdapterSeq "+params.adapter1+" "+params.adapter2+" --readFilesIn "+input.file1+" "+input.file2+" --readFilesCommand zcat --runThreadN "+str(threads)+" --outFileNamePrefix "+params.prefix+". --outSAMunmapped "+params.outsamunmapped+" --outWigType "+params.wigtype+" --outWigStrand "+params.wigstrand+" --sjdbFileChrStartEnd "+str(input.tab)+" --sjdbGTFfile "+params.gtffile+" --limitSjdbInsertNsj "+str(params.nbjuncs)+" --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate"
     # --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40"
     shell(cmd)

rule rsem:
  input: file1= "{name}.p2.Aligned.toTranscriptome.out.bam"
  output: out1="{name}.rsem.genes.results",out2="{name}.rsem.isoforms.results"
  params: rname='pl:rsem',prefix="{name}.rsem",batch='--cpus-per-task=16 --mem=32g --time=24:00:00',rsemref=config['references'][pfamily]['RSEMREF'],rsem=config['bin'][pfamily]['RSEM']
  shell: "{params.rsem}/rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  --bam --paired-end -p 16  {input.file1} {params.rsemref} {params.prefix} --time --temporary-folder /lscratch/$SLURM_JOBID --keep-intermediate-files"

rule rsemcounts:
   input: files=expand("{name}.rsem.genes.results", name=samples)
   output: "RawCountFile_RSEM_genes_filtered.txt"
   params: rname='pl:genecounts',batch='--mem=8g --time=10:00:00',dir=config['project']['workpath'],mincount=config['project']['MINCOUNTGENES'],minsamples=config['project']['MINSAMPLES'],annotate=config['references'][pfamily]['ANNOTATE']
   shell: "module load R/3.4.0_gcc-6.2.0; Rscript Scripts/rsemcounts.R '{params.dir}' '{input.files}' '{params.mincount}' '{params.minsamples}' '{params.annotate}'"

rule picard:
  input: file1= "{name}.p2.Aligned.sortedByCoord.out.bam"
  output: outstar1=temp("{name}.star_rg_added.sorted.bam"), outstar2=temp("{name}.star_rg_added.sorted.dmark.bam"),outstar3="{name}.star.duplic"
  params: rname='pl:picard',batch='--mem=24g --time=10:00:00 --gres=lscratch:800'#,picardjarpath=config['bin'][pfamily]['PICARDJARPATH']
  shell: "module load picard/1.119; java -Xmx10g  -jar $PICARDJARPATH/AddOrReplaceReadGroups.jar INPUT={input.file1} OUTPUT={output.outstar1} TMP_DIR=/lscratch/$SLURM_JOBID RGID=id RGLB=library RGPL=illumina RGPU=machine RGSM=sample; java -Xmx10g -jar $PICARDJARPATH/MarkDuplicates.jar INPUT={output.outstar1} OUTPUT={output.outstar2} TMP_DIR=/lscratch/$SLURM_JOBID CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE={output.outstar3}"

rule stats:
  input: file1= "{name}.star_rg_added.sorted.dmark.bam"
  output: outstar1="{name}.RnaSeqMetrics.txt", outstar2="{name}.flagstat.concord.txt", outstar3="{name}.InsertSizeMetrics.txt", outstar4="{name}.InsertSizeHisto.pdf"
  params: rname='pl:stats',batch='--mem=24g --time=10:00:00 --gres=lscratch:800',refflat=config['references'][pfamily]['REFFLAT'],rrnalist=config['references'][pfamily]['RRNALIST'],picardstrand=config['bin'][pfamily]['PICARDSTRAND']
  shell: "module load R/3.4.0_gcc-6.2.0;module load picard/1.119; java -Xmx10g -jar $PICARDJARPATH/CollectRnaSeqMetrics.jar REF_FLAT={params.refflat} INPUT={input.file1} OUTPUT={output.outstar1} RIBOSOMAL_INTERVALS={params.rrnalist}  STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND TMP_DIR=/lscratch/$SLURM_JOBID  VALIDATION_STRINGENCY=SILENT; java -Xmx10g -jar $PICARDJARPATH/CollectInsertSizeMetrics.jar INPUT={input.file1} OUTPUT={output.outstar3} HISTOGRAM_FILE={output.outstar4} MINIMUM_PCT=0.5 TMP_DIR=/lscratch/$SLURM_JOBID ;module load samtools; samtools flagstat {input.file1} > {output.outstar2}; module load python; python Scripts/bam_count_concord_stats.py {input.file1} >> {output.outstar2} "

rule prernaseqc:
  input: expand("{name}.star_rg_added.sorted.dmark.bam", name=samples)
  output: out1="files_to_rnaseqc.txt"
  priority: 2
  params: rname='pl:prernaseqc',batch='--mem=4g --time=04:00:00'
  run:
        with open(output.out1, "w") as out:
            out.write("Sample ID\tBam file\tNotes\n")
            for f in input:
                out.write("%s\t"  % f)
                out.write("%s\t"  % f)
                out.write("%s\n"  % f)
            out.close()

rule rnaseqc:
  input: "files_to_rnaseqc.txt"
  output: "STAR_QC"
  priority: 2
  params: rname='pl:rnaseqc',batch='--mem=24g --time=48:00:00',rrnalist=config['references'][pfamily]['RRNALIST'],gtffile=config['references'][pfamily]['GTFFILE'],genomefile=config['references'][pfamily]['GENOMEFILE']
  shell: """
         module load bwa/0.7.15
         var="{params.rrnalist}"
         if [  $var == "-" ]; then
                java -jar /usr/local/apps/rnaseqc/1.1.8/RNA-SeQC_v1.1.8.jar -n 1000 -s {input} -t {params.gtffile} -r {params.genomefile}  -o {output}
         else
                java -jar /usr/local/apps/rnaseqc/1.1.8/RNA-SeQC_v1.1.8.jar -n 1000 -s {input} -t {params.gtffile} -r {params.genomefile} -rRNA {params.rrnalist}  -o {output}
         fi
         """

rule rseqc:
    input: file1="{name}.star_rg_added.sorted.dmark.bam"
    # output: out1="{name}.strand.info",out2="{name}.inner_distance_plot.pdf",out3="{name}.GC_plot.pdf",out4="{name}.Rdist.info"
    output: out1="{name}.strand.info",out4="{name}.Rdist.info"
    # output: out1="{name}.strand.info",out2="{name}.inner_distance_plot.pdf",out3="{name}.GC_plot.pdf"
    params: bedref=config['references'][pfamily]['BEDREF'],prefix="{name}",rname="pl:rseqc"
    shell: """
           module load rseqc
    #       inner_distance.py -i {input.file1} -r {params.bedref} -o {params.prefix}
           infer_experiment.py -r {params.bedref}  -i {input.file1} > {output.out1}
    #       read_GC.py -i {input.file1}  -o {params.prefix}
           read_distribution.py -i {input.file1} -r {params.bedref} > {output.out4}
           """

#rule subread:
#   input:  file1="{name}.star_rg_added.sorted.dmark.bam",
#   output: out="{name}.star.count.info.txt", res="{name}.star.count.txt"
#   params: rname='pl:subread',batch='--time=4:00:00 --gres=lscratch:800',subreadver=config['bin'][pfamily]['SUBREADVER'],gtffile=config['references'][pfamily]['GTFFILE']
#   shell: "module load {params.subreadver}; featureCounts -T 16 -p -t exon -R -g gene_id -a {params.gtffile} --tmpDir /lscratch/$SLURM_JOBID  -o {output.out}  {input.file1}; sed '1d' {output.out} | cut -f1,7 > {output.res}"

#rule subreadoverlap:
#   input:  file1="{name}.star_rg_added.sorted.dmark.bam",
#   output: out="{name}.star.count.info.overlap.txt", res="{name}.star.count.overlap.txt"
#   params: rname='pl:subreadoverlap',batch='--cpus-per-task=16 --mem=24g --time=48:00:00 --gres=lscratch:800',subreadver=config['bin'][pfamily]['SUBREADVER'],gtffile=config['references'][pfamily]['GTFFILE']
#   shell: "module load {params.subreadver}; featureCounts -T 16 -p -t exon -R -O -g gene_id -a {params.gtffile} --tmpDir /lscratch/$SLURM_JOBID  -o {output.out}  {input.file1}; sed '1d' {output.out} | cut -f1,7 > {output.res}"

#rule genecounts: 
#   input: files=expand("{name}.star.count.txt", name=samples)
#   output: "RawCountFile_genes_filtered.txt"
#   params: rname='pl:genecounts',batch='--mem=8g --time=10:00:00',dir=config['project']['workpath'],mincount=config['project']['MINCOUNTGENES'],minsamples=config['project']['MINSAMPLES'],annotate=config['references'][pfamily]['ANNOTATE']
#   shell: "module load R/3.4.0_gcc-6.2.0; Rscript Scripts/genecounts.R '{params.dir}' '{input.files}' '{params.mincount}' '{params.minsamples}' '{params.annotate}'"

#rule junctioncounts: 
#   input: files=expand("{name}.p2.SJ.out.tab", name=samples)
#   output: "RawCountFile_junctions_filtered.txt"
#   params: rname='pl:junctioncounts',batch='--mem=8g --time=10:00:00',dir=config['project']['workpath'],mincount=config['project']['MINCOUNTJUNCTIONS'],minsamples=config['project']['MINSAMPLES']
#   shell: "module load R/3.4.0_gcc-6.2.0; Rscript Scripts/junctioncounts.R '{params.dir}' '{input.files}' '{params.mincount}' '{params.minsamples}'"

#rule genejunctioncounts: 
#   input: files=expand("{name}.p2.SJ.out.tab", name=samples)
#   output: "RawCountFile_Subread_genejunctions_filtered.txt"
#   params: rname='pl:genejunctions',batch='--mem=8g --time=10:00:00',dir=config['project']['workpath'],geneinfo=config['references'][pfamily]['GENEINFO'],mincount=config['project']['MINCOUNTGENEJUNCTIONS'],minsamples=config['project']['MINSAMPLES']
#   shell: "module load R/3.4.0_gcc-6.2.0; Rscript Scripts/genejunctioncounts.R '{params.dir}' '{input.files}' '{params.geneinfo}' '{params.mincount}' '{params.minsamples}'"

#rule joincounts:
#   input: files=expand("{name}.star.count.overlap.txt", name=samples),files2=expand("{name}.p2.ReadsPerGene.out.tab", name=samples)
#   output: "RawCountFileOverlap.txt","RawCountFileStar.txt"
#   params: rname='pl:junctioncounts',batch='--mem=8g --time=10:00:00',dir=config['project']['workpath'],starstrandcol=config['bin'][pfamily]['STARSTRANDCOL']
#   shell: "module load R/3.4.0_gcc-6.2.0; Rscript Scripts/joincounts.R '{params.dir}' '{input.files}' '{input.files2}' '{params.starstrandcol}'"