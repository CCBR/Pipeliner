from snakemake.utils import R

configfile: "run.json"

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


if pe=="yes":

   rule all:
      params: batch='--time=168:00:00'
      input: config['project']['id']+"_"+config['project']['flowcellid']+".xlsx","Reports/multiqc_report.html",
            expand("{name}.RnaSeqMetrics.txt",name=samples),"rawQC",
            "QC",expand("FQscreen/{name}_R1_001_trim_paired_screen.txt",name=samples), expand("FQscreen/{name}_R1_001_trim_paired_screen.png",name=samples), expand("FQscreen/{name}_R2_001_trim_paired_screen.txt",name=samples), expand("FQscreen/{name}_R2_001_trim_paired_screen.png",name=samples)
   rule rawfastqc:
      input: expand("{name}.R1.fastq.gz", name=samples), expand("{name}.R2.fastq.gz", name=samples)
      output: "rawQC"
      priority: 2
      params: rname='pl:rawfastqc',batch='--cpus-per-task=32 --mem=110g --time=48:00:00',fastqcver=config['bin'][pfamily]['FASTQCVER'],workpath=config['project']['workpath']
      threads: 32
      shell: "mkdir -p {output};module load {params.fastqcver}; fastqc {input} -t {threads} -o {output}; module load python/3.5; python Scripts/get_read_length.py {params.workpath}/rawQC > {output}/readlength.txt 2> {output}/readlength.err"
      
   rule trimmomatic_pe:
      input: file1= config['project']['workpath']+"/{name}.R1."+config['project']['filetype'],file2=config['project']['workpath']+"/{name}.R2."+config['project']['filetype']
      output: out11="trim/{name}_R1_001_trim_paired.fastq.gz",out12="trim/{name}_R1_001_trim_unpaired.fastq.gz",out21="trim/{name}_R2_001_trim_paired.fastq.gz",out22="trim/{name}_R2_001_trim_unpaired.fastq.gz",err="QC/{name}_run_trimmomatic.err"
      params: rname='pl:trimmomatic_pe',batch='--cpus-per-task=32 --mem=110g --time=48:00:00',trimmomaticver=config['bin'][pfamily]['TRIMMOMATICVER'],fastawithadaptersetc=config['references'][pfamily]['FASTAWITHADAPTERSETC'],seedmismatches=config['bin'][pfamily]['SEEDMISMATCHES'],palindromeclipthreshold=config['bin'][pfamily]['PALINDROMECLIPTHRESHOLD'],simpleclipthreshold=config['bin'][pfamily]['SIMPLECLIPTHRESHOLD'],leadingquality=config['bin'][pfamily]['LEADINGQUALITY'],trailingquality=config['bin'][pfamily]['TRAILINGQUALITY'],windowsize=config['bin'][pfamily]['WINDOWSIZE'],windowquality=config['bin'][pfamily]['WINDOWQUALITY'],targetlength=config['bin'][pfamily]['TARGETLENGTH'],strictness=config['bin'][pfamily]['STRICTNESS'],minlen=config['bin'][pfamily]['MINLEN'],headcroplength=config['bin'][pfamily]['HEADCROPLENGTH']
      threads:32
      shell:"module load {params.trimmomaticver}; java -classpath $TRIMMOJAR   org.usadellab.trimmomatic.TrimmomaticPE -threads {threads} {input.file1} {input.file2} {output.out11} {output.out12} {output.out21} {output.out22} ILLUMINACLIP:{params.fastawithadaptersetc}:{params.seedmismatches}:{params.palindromeclipthreshold}:{params.simpleclipthreshold} LEADING:10 TRAILING:10 MAXINFO:50:0.97 MINLEN:{params.minlen} 2> {output.err}"

   rule fastqc:
      input: expand("trim/{name}_R1_001_trim_paired.fastq.gz", name=samples), expand("trim/{name}_R2_001_trim_paired.fastq.gz", name=samples)
      output: "QC"
      priority: 2
      params: rname='pl:fastqc',batch='--cpus-per-task=32 --mem=110g --time=48:00:00',fastqcver=config['bin'][pfamily]['FASTQCVER'],workpath=config['project']['workpath']
      threads: 32
      shell: "mkdir -p {output};module load {params.fastqcver}; fastqc {input} -t {threads} -o {output}; module load python/3.5; python Scripts/get_read_length.py {params.workpath}/QC > {output}/readlength.txt  2> {output}/readlength.err"

   rule fastq_screen:
      input: file1="trim/{name}_R1_001_trim_paired.fastq.gz", file2="trim/{name}_R2_001_trim_paired.fastq.gz"
      output: out1="FQscreen/{name}_R1_001_trim_paired_screen.txt", out2="FQscreen/{name}_R1_001_trim_paired_screen.png", out3="FQscreen/{name}_R2_001_trim_paired_screen.txt", out4="FQscreen/{name}_R2_001_trim_paired_screen.png"
      params: rname='pl:fqscreen',batch='--cpus-per-task=24 --mem=64g --time=10:00:00',fastq_screen=config['bin'][pfamily]['FASTQ_SCREEN'],outdir = "FQscreen",config=config['references'][pfamily]['FASTQ_SCREEN_CONFIG']
      threads: 24
      shell:  "module load bowtie/2-2.2.9;module load perl/5.18.2; {params.fastq_screen} --conf {params.config} --outdir {params.outdir} --subset 1000000 --aligner bowtie2 --force {input.file1} {input.file2}"


   rule star1p:
      input: file1= "trim/{name}_R1_001_trim_paired.fastq.gz",file2="trim/{name}_R2_001_trim_paired.fastq.gz",qcdir="QC"#,file3="fastqc_status_checked.txt"
      output: out1= "{name}.SJ.out.tab", out3= temp("{name}.Aligned.out.bam") #,out2= "{name}.SJ.out.tab.Pass1.sjdb"
      params: rname='pl:star1p',prefix="{name}",batch='--cpus-per-task=32 --mem=110g --time=48:00:00',starver=config['bin'][pfamily]['STARVER'],stardir=config['references']['rnaseq']['STARDIR'],filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],filtertype=config['bin'][pfamily]['FILTERTYPE'],filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],adapter1=config['bin'][pfamily]['ADAPTER1'],adapter2=config['bin'][pfamily]['ADAPTER2'],
      threads: 32
      run:
        import glob,json
        rl=int(open(input.qcdir+"/readlength.txt").readlines()[0].strip())-1
        dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.stardir+'*/',recursive=False))))
        bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
        cmd="module load "+params.starver+"; STAR --genomeDir "+params.stardir+str(bestdbrl)+" --outFilterIntronMotifs "+params.filterintronmotifs+" --outSAMstrandField "+params.samstrandfield+"  --outFilterType "+params.filtertype+" --outFilterMultimapNmax "+str(params.filtermultimapnmax)+" --alignSJoverhangMin "+str(params.alignsjoverhangmin)+" --alignSJDBoverhangMin "+str(params.alignsjdboverhangmin)+"  --outFilterMismatchNmax "+str(params.filtermismatchnmax)+" --outFilterMismatchNoverLmax "+str(params.filtermismatchnoverlmax)+"  --alignIntronMin "+str(params.alignintronmin)+" --alignIntronMax "+str(params.alignintronmax)+" --alignMatesGapMax "+str(params.alignmatesgapmax)+" clip3pAdapterSeq "+params.adapter1+" "+params.adapter2+" --readFilesIn "+input.file1+" "+input.file2+" --readFilesCommand zcat --runThreadN "+str(threads)+" --outFileNamePrefix "+params.prefix+". --outSAMtype BAM Unsorted"
        shell(cmd)
        A=open("run.json",'r')
        a=eval(A.read())
        A.close()
        config=dict(a.items())
        config['project']['SJDBOVERHANG']=str(bestdbrl)
        config['project']["STARDIR"]= config['references'][pfamily]['STARDIR']+str(bestdbrl)
        with open('run.json','w') as F:
          json.dump(config, F, sort_keys = True, indent = 4,ensure_ascii=False)
        F.close()

   rule sjdb:
      input: files=expand("{name}.SJ.out.tab", name=samples)
      output: out1="uniq.filtered.SJ.out.tab"
      params: rname='pl:sjdb'
      shell: "cat {input.files} |sort|uniq|awk -F \"\\t\" '{{if ($5>0 && $6==1) {{print}}}}'|cut -f1-4|sort|uniq|grep \"^chr\"|grep -v \"^chrM\" > {output.out1}"

   rule star2p:
      input: file1= "trim/{name}_R1_001_trim_paired.fastq.gz",file2="trim/{name}_R2_001_trim_paired.fastq.gz",tab="uniq.filtered.SJ.out.tab",qcdir="QC"#, dir="STARINDEX"
      output: out1=temp("{name}.p2.Aligned.sortedByCoord.out.bam"),out2="{name}.p2.ReadsPerGene.out.tab",out3="{name}.p2.Aligned.toTranscriptome.out.bam",out4="{name}.p2.SJ.out.tab",out5="{name}.p2.Log.final.out"
      params: rname='pl:star2p',prefix="{name}.p2",batch='--cpus-per-task=32 --mem=110g --time=48:00:00',starver=config['bin'][pfamily]['STARVER'],filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],filtertype=config['bin'][pfamily]['FILTERTYPE'],filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],adapter1=config['bin'][pfamily]['ADAPTER1'],adapter2=config['bin'][pfamily]['ADAPTER2'],outsamunmapped=config['bin'][pfamily]['OUTSAMUNMAPPED'],wigtype=config['bin'][pfamily]['WIGTYPE'],wigstrand=config['bin'][pfamily]['WIGSTRAND'], gtffile=config['references'][pfamily]['GTFFILE'], nbjuncs=config['bin'][pfamily]['NBJUNCS'],stardir=config['references']['rnaseq']['STARDIR']
      threads:32
      run:
        import glob
        rl=int(open(input.qcdir+"/readlength.txt").readlines()[0].strip())-1
        dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.stardir+'*/',recursive=False))))
        bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
        cmd="module load "+params.starver+"; STAR --genomeDir "+params.stardir+str(bestdbrl)+" --outFilterIntronMotifs "+params.filterintronmotifs+" --outSAMstrandField "+params.samstrandfield+"  --outFilterType "+params.filtertype+" --outFilterMultimapNmax "+str(params.filtermultimapnmax)+" --alignSJoverhangMin "+str(params.alignsjoverhangmin)+" --alignSJDBoverhangMin "+str(params.alignsjdboverhangmin)+"  --outFilterMismatchNmax "+str(params.filtermismatchnmax)+" --outFilterMismatchNoverLmax "+str(params.filtermismatchnoverlmax)+"  --alignIntronMin "+str(params.alignintronmin)+" --alignIntronMax "+str(params.alignintronmax)+" --alignMatesGapMax "+str(params.alignmatesgapmax)+"  --clip3pAdapterSeq "+params.adapter1+" "+params.adapter2+" --readFilesIn "+input.file1+" "+input.file2+" --readFilesCommand zcat --runThreadN "+str(threads)+" --outFileNamePrefix "+params.prefix+". --outSAMunmapped "+params.outsamunmapped+" --outWigType "+params.wigtype+" --outWigStrand "+params.wigstrand+" --sjdbFileChrStartEnd "+str(input.tab)+" --sjdbGTFfile "+params.gtffile+" --limitSjdbInsertNsj "+str(params.nbjuncs)+" --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate"
        shell(cmd)

if se=="yes":

   rule all:
      params: batch='--time=168:00:00'
      input: config['project']['id']+"_"+config['project']['flowcellid']+".xlsx","Reports/multiqc_report.html",expand("{name}.RnaSeqMetrics.txt",name=samples),"rawQC","QC",expand("FQscreen/{name}_R1_001_trim_paired_screen.txt",name=samples), expand("FQscreen/{name}_R1_001_trim_paired_screen.png",name=samples)

   rule rawfastqc:
      input: expand("{name}.R1.fastq.gz", name=samples)
      output: "rawQC"
      priority: 2
      params: rname='pl:rawfastqc',batch='--cpus-per-task=32 --mem=110g --time=48:00:00',fastqcver=config['bin'][pfamily]['FASTQCVER'],workpath=config['project']['workpath']
      threads: 32
      shell: "mkdir -p {output};module load {params.fastqcver}; fastqc {input} -t {threads} -o {output}; module load python/3.5; python Scripts/get_read_length.py {params.workpath}/rawQC > {output}/readlength.txt 2> {output}/readlength.err"

   rule trimmomatic_se:
      input: file1= config['project']['workpath']+"/{name}.R1."+config['project']['filetype']
      output: out11=temp("trim/{name}_R1_001_trim_paired.fastq.gz"),err="QC/{name}_run_trimmomatic.err"
      params: rname='pl:trimmomatic_se',batch='--cpus-per-task=32 --mem=110g --time=48:00:00',trimmomaticver=config['bin'][pfamily]['TRIMMOMATICVER'],fastawithadaptersetc=config['references'][pfamily]['FASTAWITHADAPTERSETC'],seedmismatches=config['bin'][pfamily]['SEEDMISMATCHES'],palindromeclipthreshold=config['bin'][pfamily]['PALINDROMECLIPTHRESHOLD'],simpleclipthreshold=config['bin'][pfamily]['SIMPLECLIPTHRESHOLD'],leadingquality=config['bin'][pfamily]['LEADINGQUALITY'],trailingquality=config['bin'][pfamily]['TRAILINGQUALITY'],windowsize=config['bin'][pfamily]['WINDOWSIZE'],windowquality=config['bin'][pfamily]['WINDOWQUALITY'],targetlength=config['bin'][pfamily]['TARGETLENGTH'],strictness=config['bin'][pfamily]['STRICTNESS'],minlen=config['bin'][pfamily]['MINLEN'],headcroplength=config['bin'][pfamily]['HEADCROPLENGTH']
      threads:32
      shell:"module load {params.trimmomaticver}; java -classpath $TRIMMOJAR   org.usadellab.trimmomatic.TrimmomaticSE -threads {threads} {input.file1} {output.out11} ILLUMINACLIP:{params.fastawithadaptersetc}:{params.seedmismatches}:{params.palindromeclipthreshold}:{params.simpleclipthreshold} LEADING:10 TRAILING:10 MAXINFO:50:0.97 MINLEN:{params.minlen} 2> {output.err}"

   rule fastqc:
      input: expand("trim/{name}_R1_001_trim_paired.fastq.gz", name=samples)
      output: "QC"
      priority: 2
      params: rname='pl:fastqc',batch='--cpus-per-task=32 --mem=110g --time=48:00:00',fastqcver=config['bin'][pfamily]['FASTQCVER'],workpath=config['project']['workpath']
      threads: 32
      shell: "mkdir -p {output};module load {params.fastqcver}; fastqc {input} -t {threads} -o {output}; module load python/3.5; python Scripts/get_read_length.py {params.workpath}/QC > {output}/readlength.txt  2> {output}/readlength.err"

   rule fastq_screen:
      input: file1="trim/{name}_R1_001_trim_paired.fastq.gz"
      output: out1="FQscreen/{name}_R1_001_trim_paired_screen.txt", out2="FQscreen/{name}_R1_001_trim_paired_screen.png"
      params: rname='pl:fqscreen',batch='--cpus-per-task=24 --mem=64g --time=10:00:00',fastq_screen=config['bin'][pfamily]['FASTQ_SCREEN'],outdir = "FQscreen",config=config['references'][pfamily]['FASTQ_SCREEN_CONFIG']
      threads: 24
      shell:  "module load bowtie/2-2.2.9;module load perl/5.18.2; {params.fastq_screen} --conf {params.config} --outdir {params.outdir} --subset 1000000 --aligner bowtie2 --force {input.file1}"

   rule star1p:
      input: file1= "trim/{name}_R1_001_trim_paired.fastq.gz",qcdir="QC"#,file3="fastqc_status_checked.txt"
      output: out1= "{name}.SJ.out.tab", out3= temp("{name}.Aligned.out.bam") #,out2= "{name}.SJ.out.tab.Pass1.sjdb"
      params: rname='pl:star1p',prefix="{name}",batch='--cpus-per-task=32 --mem=110g --time=48:00:00',starver=config['bin'][pfamily]['STARVER'],stardir=config['references']['rnaseq']['STARDIR'],filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],filtertype=config['bin'][pfamily]['FILTERTYPE'],filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],adapter1=config['bin'][pfamily]['ADAPTER1'],adapter2=config['bin'][pfamily]['ADAPTER2'],
      threads: 32
      run:
        import glob,json
        rl=int(open(input.qcdir+"/readlength.txt").readlines()[0].strip())-1
        dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.stardir+'*/',recursive=False))))
        bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
        cmd="module load "+params.starver+"; STAR --genomeDir "+params.stardir+str(bestdbrl)+" --outFilterIntronMotifs "+params.filterintronmotifs+" --outSAMstrandField "+params.samstrandfield+"  --outFilterType "+params.filtertype+" --outFilterMultimapNmax "+str(params.filtermultimapnmax)+" --alignSJoverhangMin "+str(params.alignsjoverhangmin)+" --alignSJDBoverhangMin "+str(params.alignsjdboverhangmin)+"  --outFilterMismatchNmax "+str(params.filtermismatchnmax)+" --outFilterMismatchNoverLmax "+str(params.filtermismatchnoverlmax)+"  --alignIntronMin "+str(params.alignintronmin)+" --alignIntronMax "+str(params.alignintronmax)+" --alignMatesGapMax "+str(params.alignmatesgapmax)+" clip3pAdapterSeq "+params.adapter1+" "+params.adapter2+" --readFilesIn "+input.file1+" --readFilesCommand zcat --runThreadN "+str(threads)+" --outFileNamePrefix "+params.prefix+". --outSAMtype BAM Unsorted"
        shell(cmd)
        A=open("run.json",'r')
        a=eval(A.read())
        A.close()
        config=dict(a.items())
        config['project']['SJDBOVERHANG']=str(bestdbrl)
        config['project']["STARDIR"]= config['references'][pfamily]['STARDIR']+str(bestdbrl)
        with open('run.json','w') as F:
          json.dump(config, F, sort_keys = True, indent = 4,ensure_ascii=False)
        F.close()

   rule sjdb:
      input: files=expand("{name}.SJ.out.tab", name=samples)
      output: out1="uniq.filtered.SJ.out.tab"
      params: rname='pl:sjdb'
      shell: "cat {input.files} |sort|uniq|awk -F \"\\t\" '{{if ($5>0 && $6==1) {{print}}}}'|cut -f1-4|sort|uniq|grep \"^chr\"|grep -v \"^chrM\" > {output.out1}"

   rule star2p:
      input: file1= "trim/{name}_R1_001_trim_paired.fastq.gz",tab="uniq.filtered.SJ.out.tab",qcdir="QC"#, dir="STARINDEX"
      output: out1=temp("{name}.p2.Aligned.sortedByCoord.out.bam"),out2="{name}.p2.ReadsPerGene.out.tab",out3="{name}.p2.Aligned.toTranscriptome.out.bam",out4="{name}.p2.SJ.out.tab",out5="{name}.p2.Log.final.out"
      params: rname='pl:star2p',prefix="{name}.p2",batch='--cpus-per-task=32 --mem=110g --time=48:00:00',starver=config['bin'][pfamily]['STARVER'],filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],filtertype=config['bin'][pfamily]['FILTERTYPE'],filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],adapter1=config['bin'][pfamily]['ADAPTER1'],adapter2=config['bin'][pfamily]['ADAPTER2'],outsamunmapped=config['bin'][pfamily]['OUTSAMUNMAPPED'],wigtype=config['bin'][pfamily]['WIGTYPE'],wigstrand=config['bin'][pfamily]['WIGSTRAND'], gtffile=config['references'][pfamily]['GTFFILE'], nbjuncs=config['bin'][pfamily]['NBJUNCS'],stardir=config['references']['rnaseq']['STARDIR']
      threads:32
      run:
        import glob
        rl=int(open(input.qcdir+"/readlength.txt").readlines()[0].strip())-1
        dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.stardir+'*/',recursive=False))))
        bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
        cmd="module load "+params.starver+"; STAR --genomeDir "+params.stardir+str(bestdbrl)+" --outFilterIntronMotifs "+params.filterintronmotifs+" --outSAMstrandField "+params.samstrandfield+"  --outFilterType "+params.filtertype+" --outFilterMultimapNmax "+str(params.filtermultimapnmax)+" --alignSJoverhangMin "+str(params.alignsjoverhangmin)+" --alignSJDBoverhangMin "+str(params.alignsjdboverhangmin)+"  --outFilterMismatchNmax "+str(params.filtermismatchnmax)+" --outFilterMismatchNoverLmax "+str(params.filtermismatchnoverlmax)+"  --alignIntronMin "+str(params.alignintronmin)+" --alignIntronMax "+str(params.alignintronmax)+" --alignMatesGapMax "+str(params.alignmatesgapmax)+"  --clip3pAdapterSeq "+params.adapter1+" "+params.adapter2+" --readFilesIn "+input.file1+" --readFilesCommand zcat --runThreadN "+str(threads)+" --outFileNamePrefix "+params.prefix+". --outSAMunmapped "+params.outsamunmapped+" --outWigType "+params.wigtype+" --outWigStrand "+params.wigstrand+" --sjdbFileChrStartEnd "+str(input.tab)+" --sjdbGTFfile "+params.gtffile+" --limitSjdbInsertNsj "+str(params.nbjuncs)+" --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate"
        shell(cmd)



rule picard:
   input: file1= "{name}.p2.Aligned.sortedByCoord.out.bam"
   output: outstar1=temp("{name}.star_rg_added.sorted.bam"), outstar2="{name}.star_rg_added.sorted.dmark.bam",outstar3="{name}.star.duplic"
   params: rname='pl:picard',batch='--mem=24g --time=10:00:00 --gres=lscratch:800',picardver=config['bin'][pfamily]['PICARDVER']#,picardjarpath=config['bin'][pfamily]['PICARDJARPATH']
   shell: "module load {params.picardver}; java -Xmx10g  -jar $PICARDJARPATH/AddOrReplaceReadGroups.jar INPUT={input.file1} OUTPUT={output.outstar1} TMP_DIR=/lscratch/$SLURM_JOBID RGID=id RGLB=library RGPL=illumina RGPU=machine RGSM=sample; java -Xmx10g -jar $PICARDJARPATH/MarkDuplicates.jar INPUT={output.outstar1} OUTPUT={output.outstar2} TMP_DIR=/lscratch/$SLURM_JOBID CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE={output.outstar3}"

rule stats:
   input: file1= "{name}.star_rg_added.sorted.dmark.bam"
   output: outstar1="{name}.RnaSeqMetrics.txt", outstar2="{name}.flagstat.concord.txt"
   params: rname='pl:stats',batch='--mem=24g --time=10:00:00 --gres=lscratch:800',picardver=config['bin'][pfamily]['PICARDVER'],refflat=config['references'][pfamily]['REFFLAT'],rrnalist=config['references'][pfamily]['RRNALIST'],picardstrand=config['bin'][pfamily]['PICARDSTRAND']
   shell: "module load R/3.4.0_gcc-6.2.0;module load {params.picardver}; java -Xmx10g -jar $PICARDJARPATH/CollectRnaSeqMetrics.jar REF_FLAT={params.refflat} INPUT={input.file1} OUTPUT={output.outstar1} RIBOSOMAL_INTERVALS={params.rrnalist}  STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND TMP_DIR=/lscratch/$SLURM_JOBID  VALIDATION_STRINGENCY=SILENT;module load samtools; samtools flagstat {input.file1} > {output.outstar2}; module load python; python Scripts/bam_count_concord_stats.py {input.file1} >> {output.outstar2} "


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
   input: expand("{name}.Rdist.info",name=samples),expand("FQscreen/{name}_R1_001_trim_paired_screen.png",name=samples),expand("{name}.flagstat.concord.txt",name=samples),expand("{name}.RnaSeqMetrics.txt",name=samples)
   output: "Reports/multiqc_report.html"
   params: rname="pl:multiqc",multiqc=config['bin'][pfamily]['MULTIQC'],qcconfig=config['bin'][pfamily]['CONFMULTIQC']
   threads: 1
   shell:  """
            module load {params.multiqc}
            cd Reports && multiqc -f -c {params.qcconfig}  ../
            """

if pe=="yes":

   rule RNAseq_generate_QC_table:
      input: expand("QC/{name}_run_trimmomatic.err",name=samples), expand("{name}.star.duplic",name=samples), expand("{name}.p2.Log.final.out",name=samples), expand("{name}.RnaSeqMetrics.txt",name=samples)
      output: config['project']['id']+"_"+config['project']['flowcellid']+".xlsx"
      params: project=config['project']['id'],flowcell=config['project']['flowcellid'],dir=config['project']['workpath'],rname="pl:QC_table"
      shell: "module load perl/5.18.2; perl Scripts/CollectPipelineStats2Tab_v2.3.pl -p {params.project} -f {params.flowcell} -d {params.dir} -r 5 -e 2; perl Scripts/Tab2Excel_v2.3.pl -i {params.project}_{params.flowcell} -r 5"

if se=="yes":

   rule RNAseq_generate_QC_table:
      input: expand("QC/{name}_run_trimmomatic.err",name=samples), expand("{name}.star.duplic",name=samples), expand("{name}.p2.Log.final.out",name=samples), expand("{name}.RnaSeqMetrics.txt",name=samples)
      output: config['project']['id']+"_"+config['project']['flowcellid']+".xlsx"
      params: project=config['project']['id'],flowcell=config['project']['flowcellid'],dir=config['project']['workpath'],rname="pl:QC_table"
      shell: "module load perl/5.18.2; perl Scripts/CollectPipelineStats2Tab_v2.3.pl -p {params.project} -f {params.flowcell} -d {params.dir} -r 5 -e 1; perl Scripts/Tab2Excel_v2.3.pl -i {params.project}_{params.flowcell} -r 5"



rule rseqc:
   input: file1="{name}.star_rg_added.sorted.dmark.bam"
   output: out1="{name}.strand.info",out4="{name}.Rdist.info"
   params: bedref=config['references'][pfamily]['BEDREF'],prefix="{name}",rname="pl:rseqc"
   shell: """
           module load rseqc
           infer_experiment.py -r {params.bedref}  -i {input.file1} > {output.out1}
           read_distribution.py -i {input.file1} -r {params.bedref} > {output.out4}
           """
