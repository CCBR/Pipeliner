from snakemake.utils import R

configfile: "run.json"

if config['project']['DEG'] == "yes" and config['project']['TRIM'] == "yes":
  rule all:
     params: batch='--time=168:00:00'
     input: config['project']['id']+"_"+config['project']['flowcellid']+".xlsx","Reports/multiqc_report.html",
            expand("{name}.RnaSeqMetrics.txt",name=samples),"rawQC",
            "QC",expand("FQscreen/{name}_R1_001_trim_paired_screen.txt",name=samples), expand("FQscreen/{name}_R1_001_trim_paired_screen.png",name=samples), expand("FQscreen/{name}_R2_001_trim_paired_screen.txt",name=samples), expand("FQscreen/{name}_R2_001_trim_paired_screen.png",name=samples)


elif config['project']['DEG'] == "no" and config['project']['TRIM'] == "yes":
  rule all:
     params: batch='--time=168:00:00'
     input: config['project']['id']+"_"+config['project']['flowcellid']+".xlsx","Reports/multiqc_report.html",
            expand("{name}.RnaSeqMetrics.txt",name=samples),"rawQC",
            "QC",expand("FQscreen/{name}_R1_001_trim_paired_screen.txt",name=samples), expand("FQscreen/{name}_R1_001_trim_paired_screen.png",name=samples), expand("FQscreen/{name}_R2_001_trim_paired_screen.txt",name=samples), expand("FQscreen/{name}_R2_001_trim_paired_screen.png",name=samples)


elif config['project']['DEG'] == "yes" and config['project']['TRIM'] == "no":
  rule all:
     input: config['project']['id']+"_"+config['project']['flowcellid']+".xlsx","Reports/multiqc_report.html",
            expand("{name}.RnaSeqMetrics.txt",name=samples),"QC",expand("FQscreen/{name}.R1_screen.txt",name=samples),expand("FQscreen/{name}.R1_screen.png",name=samples),expand("FQscreen/{name}.R2_screen.txt",name=samples),expand("FQscreen/{name}.R2_screen.png",name=samples),expand("{name}.InsertSizeMetrics.txt",name=samples)
     params: batch='--time=168:00:00'
else:
  rule all:
     params: batch='--time=168:00:00'
     input: config['project']['id']+"_"+config['project']['flowcellid']+".xlsx","Reports/multiqc_report.html",
            expand("{name}.RnaSeqMetrics.txt",name=samples),"QC",expand("FQscreen/{name}.R1_screen.txt",name=samples),expand("FQscreen/{name}.R1_screen.png",name=samples),expand("FQscreen/{name}.R2_screen.txt",name=samples),expand("FQscreen/{name}.R2_screen.png",name=samples),expand("{name}.InsertSizeMetrics.txt",name=samples)

# rule fastq_screen:
#       input: "{name}.R1.fastq.gz","{name}.R2.fastq.gz"
#       output: "FQscreen/{name}.R1_screen.txt","FQscreen/{name}.R1_screen.png","FQscreen/{name}.R2_screen.txt","FQscreen/{name}.R2_screen.png"
#       params: rname='pl:fqscreen',fastq_screen=config['bin'][pfamily]['FASTQ_SCREEN'],outdir = "FQscreen",config=config['references'][pfamily]['FASTQ_SCREEN_CONFIG']
#       threads: 24
#       shell:  "module load bowtie/2-2.2.9;module load perl/5.18.2; {params.fastq_screen}_v2 --conf {params.config} --outdir {params.outdir} --subset 1000000 --aligner bowtie2 --force {input}"

rule rawfastqc:
      input: expand("{name}.R1.fastq.gz", name=samples), expand("{name}.R2.fastq.gz", name=samples)
      output: "rawQC"
      priority: 2
      params: rname='pl:rawfastqc',batch='--cpus-per-task=32 --mem=110g --time=48:00:00',fastqcver=config['bin'][pfamily]['FASTQCVER'],workpath=config['project']['workpath']
      threads: 32
      shell: "mkdir -p {output};module load {params.fastqcver}; fastqc {input} -t {threads} -o {output}; module load python/3.5; python Scripts/get_read_length.py {params.workpath}/rawQC > {output}/readlength.txt 2> {output}/readlength.err"

if config['project']['TRIM'] == "yes":

   rule trimmomatic_pe:
      input: file1= config['project']['workpath']+"/{name}.R1."+config['project']['filetype'],file2=config['project']['workpath']+"/{name}.R2."+config['project']['filetype']
      output: out11=temp("trim/{name}_R1_001_trim_paired.fastq.gz"),out12=temp("trim/{name}_R1_001_trim_unpaired.fastq.gz"),out21=temp("trim/{name}_R2_001_trim_paired.fastq.gz"),out22=temp("trim/{name}_R2_001_trim_unpaired.fastq.gz"),err="QC/{name}_run_trimmomatic.err"
      params: rname='pl:trimmomatic_pe',batch='--cpus-per-task=32 --mem=110g --time=48:00:00',trimmomaticver=config['bin'][pfamily]['TRIMMOMATICVER'],fastawithadaptersetc=config['references'][pfamily]['FASTAWITHADAPTERSETC'],seedmismatches=config['bin'][pfamily]['SEEDMISMATCHES'],palindromeclipthreshold=config['bin'][pfamily]['PALINDROMECLIPTHRESHOLD'],simpleclipthreshold=config['bin'][pfamily]['SIMPLECLIPTHRESHOLD'],leadingquality=config['bin'][pfamily]['LEADINGQUALITY'],trailingquality=config['bin'][pfamily]['TRAILINGQUALITY'],windowsize=config['bin'][pfamily]['WINDOWSIZE'],windowquality=config['bin'][pfamily]['WINDOWQUALITY'],targetlength=config['bin'][pfamily]['TARGETLENGTH'],strictness=config['bin'][pfamily]['STRICTNESS'],minlen=config['bin'][pfamily]['MINLEN'],headcroplength=config['bin'][pfamily]['HEADCROPLENGTH']
      threads:32
      shell:"module load {params.trimmomaticver}; java -classpath $TRIMMOJAR   org.usadellab.trimmomatic.TrimmomaticPE -threads {threads} {input.file1} {input.file2} {output.out11} {output.out12} {output.out21} {output.out22} ILLUMINACLIP:{params.fastawithadaptersetc}:{params.seedmismatches}:{params.palindromeclipthreshold}:{params.simpleclipthreshold} MINLEN:{params.minlen} 2> {output.err}"

   rule fastqc:
      input: expand("trim/{name}_R1_001_trim_paired.fastq.gz", name=samples), expand("trim/{name}_R2_001_trim_paired.fastq.gz", name=samples)
      output: "QC"
      priority: 2
      params: rname='pl:fastqc',batch='--cpus-per-task=32 --mem=110g --time=48:00:00',fastqcver=config['bin'][pfamily]['FASTQCVER'],workpath=config['project']['workpath']
      threads: 32
      shell: "mkdir -p {output};module load {params.fastqcver}; fastqc {input} -t {threads} -o {output}; module load python/3.5; python Scripts/get_read_length.py {params.workpath}/QC > {output}/readlength.txt  2> {output}/readlength.err"

   rule fastq_screen:
      input: expand("trim/{name}_R1_001_trim_paired.fastq.gz", name=samples), expand("trim/{name}_R2_001_trim_paired.fastq.gz", name=samples)
      output: expand("FQscreen/{name}_R1_001_trim_paired_screen.txt",name=samples), expand("FQscreen/{name}_R1_001_trim_paired_screen.png",name=samples), expand("FQscreen/{name}_R2_001_trim_paired_screen.txt",name=samples), expand("FQscreen/{name}_R2_001_trim_paired_screen.png",name=samples)
      params: rname='pl:fqscreen',fastq_screen=config['bin'][pfamily]['FASTQ_SCREEN'],outdir = "FQscreen",config=config['references'][pfamily]['FASTQ_SCREEN_CONFIG']
      threads: 24
      shell:  "module load bowtie/2-2.2.9;module load perl/5.18.2; {params.fastq_screen} --conf {params.config} --outdir {params.outdir} --subset 1000000 --aligner bowtie2 --force {input}"


   rule star1p:
      input: file1= "trim/{name}_R1_001_trim_paired.fastq.gz",file2="trim/{name}_R2_001_trim_paired.fastq.gz",qcdir="QC"#,file3="fastqc_status_checked.txt"
      output: out1= "{name}.SJ.out.tab", out3= temp("{name}.Aligned.out.bam") #,out2= "{name}.SJ.out.tab.Pass1.sjdb"
      params: rname='pl:star1p',prefix="{name}",batch='--cpus-per-task=32 --mem=110g --time=48:00:00',starver=config['bin'][pfamily]['STARVER'],stardir=config['references']['rnaseqvargerm']['STARREF'],filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],filtertype=config['bin'][pfamily]['FILTERTYPE'],filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],adapter1=config['bin'][pfamily]['ADAPTER1'],adapter2=config['bin'][pfamily]['ADAPTER2'],
      threads: 32
      run:
        import glob,json
        rl=int(open(input.qcdir+"/readlength.txt").readlines()[0].strip())-1
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
      input: file1= "trim/{name}_R1_001_trim_paired.fastq.gz",file2="trim/{name}_R2_001_trim_paired.fastq.gz",tab=expand("{name}.SJ.out.tab",name=samples),qcdir="QC"#, dir="STARINDEX"
      output: out1=temp("{name}.p2.Aligned.sortedByCoord.out.bam"),out2="{name}.p2.ReadsPerGene.out.tab",out3="{name}.p2.Aligned.toTranscriptome.out.bam",out4="{name}.p2.SJ.out.tab",out5="{name}.p2.Log.final.out"
      params: rname='pl:star2p',prefix="{name}.p2",batch='--cpus-per-task=32 --mem=110g --time=48:00:00',starver=config['bin'][pfamily]['STARVER'],filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],filtertype=config['bin'][pfamily]['FILTERTYPE'],filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],adapter1=config['bin'][pfamily]['ADAPTER1'],adapter2=config['bin'][pfamily]['ADAPTER2'],outsamunmapped=config['bin'][pfamily]['OUTSAMUNMAPPED'],wigtype=config['bin'][pfamily]['WIGTYPE'],wigstrand=config['bin'][pfamily]['WIGSTRAND'], gtffile=config['references'][pfamily]['GTFFILE2'], nbjuncs=config['bin'][pfamily]['NBJUNCS'],stardir=config['references']['rnaseqvargerm']['STARREF']
      threads:32
      run:
        import glob
        rl=int(open(input.qcdir+"/readlength.txt").readlines()[0].strip())-1
        dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.stardir+'*/',recursive=False))))
        bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
        cmd="module load "+params.starver+"; STAR --genomeDir "+params.stardir+str(bestdbrl)+" --outFilterIntronMotifs "+params.filterintronmotifs+" --outSAMstrandField "+params.samstrandfield+"  --outFilterType "+params.filtertype+" --outFilterMultimapNmax "+str(params.filtermultimapnmax)+" --alignSJoverhangMin "+str(params.alignsjoverhangmin)+" --alignSJDBoverhangMin "+str(params.alignsjdboverhangmin)+"  --outFilterMismatchNmax "+str(params.filtermismatchnmax)+" --outFilterMismatchNoverLmax "+str(params.filtermismatchnoverlmax)+"  --alignIntronMin "+str(params.alignintronmin)+" --alignIntronMax "+str(params.alignintronmax)+" --alignMatesGapMax "+str(params.alignmatesgapmax)+"  --clip3pAdapterSeq "+params.adapter1+" "+params.adapter2+" --readFilesIn "+input.file1+" "+input.file2+" --readFilesCommand zcat --runThreadN "+str(threads)+" --outFileNamePrefix "+params.prefix+". --outSAMunmapped "+params.outsamunmapped+" --outWigType "+params.wigtype+" --outWigStrand "+params.wigstrand+" --sjdbFileChrStartEnd "+str(input.tab)+" --sjdbGTFfile "+params.gtffile+" --limitSjdbInsertNsj "+str(params.nbjuncs)+" --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate"
         # --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 40"
        shell(cmd)


else:


   rule fastqc2:
      input: expand(config['project']['workpath']+"/{name}.R1."+config['project']['filetype'], name=samples), expand(config['project']['workpath']+"/{name}.R2."+config['project']['filetype'], name=samples)
      output: "QC"
      priority: 2
      params: rname='pl:fastqc',batch='--cpus-per-task=32 --mem=110g --time=48:00:00',fastqcver=config['bin'][pfamily]['FASTQCVER']
      threads: 32
      shell: "mkdir -p {output};module load {params.fastqcver}; fastqc {input} -t {threads} -o {output}"

   rule star1p:
      input: file1= config['project']['workpath']+"/{name}.R1."+config['project']['filetype'],file2=config['project']['workpath']+"/{name}.R2."+config['project']['filetype']
      output: out1= "{name}.SJ.out.tab"#,out2= "{name}.SJ.out.tab.Pass1.sjdb"
      params: rname='pl:star1p',prefix="{name}",batch='--cpus-per-task=32 --mem=110g --time=48:00:00',starver=config['bin'][pfamily]['STARVER'],stardir=config['references'][pfamily]['STARDIR']+config['project']["SJDBOVERHANG"],filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],filtertype=config['bin'][pfamily]['FILTERTYPE'],filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],adapter1=config['bin'][pfamily]['ADAPTER1'],adapter2=config['bin'][pfamily]['ADAPTER2']
      threads: 32
      shell:"module load {params.starver}; STAR --genomeDir {params.stardir} --outFilterIntronMotifs {params.filterintronmotifs} --outSAMstrandField {params.samstrandfield}  --outFilterType {params.filtertype} --outFilterMultimapNmax {params.filtermultimapnmax} --alignSJoverhangMin {params.alignsjoverhangmin} --alignSJDBoverhangMin {params.alignsjdboverhangmin}  --outFilterMismatchNmax {params.filtermismatchnmax} --outFilterMismatchNoverLmax {params.filtermismatchnoverlmax}  --alignIntronMin {params.alignintronmin} --alignIntronMax {params.alignintronmax} --alignMatesGapMax {params.alignmatesgapmax} clip3pAdapterSeq {params.adapter1} {params.adapter2} --readFilesIn {input.file1} {input.file2} --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {params.prefix}. --outSAMtype BAM Unsorted"

   rule star2p:
      input: file1= config['project']['workpath']+"/{name}.R1."+config['project']['filetype'],file2=config['project']['workpath']+"/{name}.R2."+config['project']['filetype'],tab=expand("{name}.SJ.out.tab",name=samples)#,dir=config['bin'][pfamily]['STARINDEX']
      output: out1="{name}.p2.Aligned.sortedByCoord.out.bam", out2="{name}.p2.ReadsPerGene.out.tab", out3="{name}.p2.Aligned.toTranscriptome.out.bam",out4="{name}.p2.SJ.out.tab" #"{name}.p2.Aligned.out.sam"
      params: rname='pl:star2p',prefix="{name}.p2",batch='--cpus-per-task=32 --mem=110g --time=48:00:00',starver=config['bin'][pfamily]['STARVER'],filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],filtertype=config['bin'][pfamily]['FILTERTYPE'],filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],adapter1=config['bin'][pfamily]['ADAPTER1'],adapter2=config['bin'][pfamily]['ADAPTER2'],outsamunmapped=config['bin'][pfamily]['OUTSAMUNMAPPED'],wigtype=config['bin'][pfamily]['WIGTYPE'],wigstrand=config['bin'][pfamily]['WIGSTRAND'], gtffile=config['references'][pfamily]['GTFFILE'], nbjuncs=config['bin'][pfamily]['NBJUNCS'],stardir=config['references'][pfamily]['STARDIR']+config['project']["SJDBOVERHANG"]
      threads:32
      shell:"module load {params.starver}; STAR --genomeDir {params.stardir} --outFilterIntronMotifs {params.filterintronmotifs} --outSAMstrandField {params.samstrandfield}  --outFilterType {params.filtertype} --outFilterMultimapNmax {params.filtermultimapnmax} --alignSJoverhangMin {params.alignsjoverhangmin} --alignSJDBoverhangMin {params.alignsjdboverhangmin}  --outFilterMismatchNmax {params.filtermismatchnmax} --outFilterMismatchNoverLmax {params.filtermismatchnoverlmax}  --alignIntronMin {params.alignintronmin} --alignIntronMax {params.alignintronmax} --alignMatesGapMax {params.alignmatesgapmax}  --clip3pAdapterSeq {params.adapter1} {params.adapter2} --readFilesIn {input.file1} {input.file2} --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {params.prefix}. --outSAMunmapped {params.outsamunmapped} --outWigType {params.wigtype} --outWigStrand {params.wigstrand} --sjdbFileChrStartEnd {input.tab} --sjdbGTFfile {params.gtffile} --limitSjdbInsertNsj {params.nbjuncs} --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate"


rule rsem:
  input: file1= "{name}.p2.Aligned.toTranscriptome.out.bam"
  output: out1="{name}.rsem.genes.results",out2="{name}.rsem.isoforms.results"
  params: rname='pl:rsem',prefix="{name}.rsem",batch='--cpus-per-task=16 --mem=32g --time=24:00:00',rsemref=config['references'][pfamily]['RSEMREF'],rsem=config['bin'][pfamily]['RSEM']
  shell: "{params.rsem}/rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  --bam --paired-end -p 16  {input.file1} {params.rsemref} {params.prefix}"

rule picard:
  input: file1= "{name}.p2.Aligned.sortedByCoord.out.bam"
  output: outstar1=temp("{name}.star_rg_added.sorted.bam"), outstar2="{name}.star_rg_added.sorted.dmark.bam",outstar3="{name}.star.duplic"
  params: rname='pl:picard',batch='--mem=24g --time=10:00:00 --gres=lscratch:800',picardver=config['bin'][pfamily]['PICARDVER']#,picardjarpath=config['bin'][pfamily]['PICARDJARPATH']
  shell: "module load {params.picardver}; java -Xmx10g  -jar $PICARDJARPATH/AddOrReplaceReadGroups.jar INPUT={input.file1} OUTPUT={output.outstar1} TMP_DIR=/lscratch/$SLURM_JOBID RGID=id RGLB=library RGPL=illumina RGPU=machine RGSM=sample; java -Xmx10g -jar $PICARDJARPATH/MarkDuplicates.jar INPUT={output.outstar1} OUTPUT={output.outstar2} TMP_DIR=/lscratch/$SLURM_JOBID CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE={output.outstar3}"

rule stats:
  input: file1= "{name}.star_rg_added.sorted.dmark.bam"
  output: outstar1="{name}.RnaSeqMetrics.txt", outstar2="{name}.flagstat.concord.txt", outstar3="{name}.InsertSizeMetrics.txt", outstar4="{name}.InsertSizeHisto.pdf"
  params: rname='pl:stats',batch='--mem=24g --time=10:00:00 --gres=lscratch:800',picardver=config['bin'][pfamily]['PICARDVER'],refflat=config['references'][pfamily]['REFFLAT2'],rrnalist=config['references'][pfamily]['RRNALIST2'],picardstrand=config['bin'][pfamily]['PICARDSTRAND']
  shell: "module load R;module load {params.picardver}; java -Xmx10g -jar $PICARDJARPATH/CollectRnaSeqMetrics.jar REF_FLAT={params.refflat} INPUT={input.file1} OUTPUT={output.outstar1} RIBOSOMAL_INTERVALS={params.rrnalist}  STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND TMP_DIR=/lscratch/$SLURM_JOBID  VALIDATION_STRINGENCY=SILENT; java -Xmx10g -jar $PICARDJARPATH/CollectInsertSizeMetrics.jar INPUT={input.file1} OUTPUT={output.outstar3} HISTOGRAM_FILE={output.outstar4} MINIMUM_PCT=0.5 TMP_DIR=/lscratch/$SLURM_JOBID ;module load samtools; samtools flagstat {input.file1} > {output.outstar2}; module load python; python Scripts/bam_count_concord_stats.py {input.file1} >> {output.outstar2} "


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
  params: rname='pl:rnaseqc',batch='--mem=24g --time=48:00:00',bwaver=config['bin'][pfamily]['BWAVER'],rrnalist=config['references'][pfamily]['RRNALIST'],rnaseqcver=config['bin'][pfamily]['RNASEQCVER'],gtffile=config['references'][pfamily]['GTFFILE'],genomefile=config['references'][pfamily]['GENOMEFILE']
  shell: """
         module load {params.bwaver}
         var="{params.rrnalist}"
         if [  $var == "-" ]; then
                java -jar {params.rnaseqcver} -n 1000 -s {input} -t {params.gtffile} -r {params.genomefile}  -o {output}
         else
                java -jar {params.rnaseqcver} -n 1000 -s {input} -t {params.gtffile} -r {params.genomefile} -rRNA {params.rrnalist}  -o {output}
         fi
         """

rule rnaseq_multiqc:
    input: expand("{name}.Rdist.info",name=samples),expand("FQscreen/{name}_R1_001_trim_paired_screen.png",name=samples),expand("FQscreen/{name}_R2_001_trim_paired_screen.png",name=samples),expand("{name}.flagstat.concord.txt",name=samples),expand("{name}.RnaSeqMetrics.txt",name=samples),expand("{name}.InsertSizeMetrics.txt",name=samples)
    output: "Reports/multiqc_report.html"
#    params: rname="pl:multiqc",pythonpath=config['bin'][pfamily]['PYTHONPATH'],multiqc=config['bin'][pfamily]['MULTIQC']
    params: rname="pl:multiqc",multiqc=config['bin'][pfamily]['MULTIQC'],qcconfig=config['bin'][pfamily]['CONFMULTIQC']
    threads: 1
    shell:  """
            module load {params.multiqc}
            # cd Reports && multiqc -f -e featureCounts -e picard ../
#            cd Reports && multiqc -f -e featureCounts ../
            cd Reports && multiqc -f -c {params.qcconfig}  ../
            """

rule RNAseq_generate_QC_table:
    input: expand("QC/{name}_run_trimmomatic.err",name=samples), expand("{name}.star.duplic",name=samples), expand("{name}.p2.Log.final.out",name=samples), expand("{name}.RnaSeqMetrics.txt",name=samples)
    output: config['project']['id']+"_"+config['project']['flowcellid']+".xlsx"
#    output: "QC_table.xlsx"
    params: project=config['project']['id'],flowcell=config['project']['flowcellid'],dir=config['project']['workpath'],rname="pl:QC_table"
    shell: "module load perl/5.18.2; perl Scripts/CollectPipelineStats2Tab_v2.3.pl -p {params.project} -f {params.flowcell} -d {params.dir} -r 5 -e 2; perl Scripts/Tab2Excel_v2.3.pl -i {params.project}_{params.flowcell} -r 5"

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
