from snakemake.utils import R

configfile: "run.json"

#samples=config['project']['units']
#degsamples=config['project']['contrasts']['rsamps']

if config['project']['DEG'] == "yes" and config['project']['TRIM'] == "yes":
  rule all:
     params: batch='--time=168:00:00'
     input: "STAR_QC","Reports/multiqc_report.html",
            "ebseq_completed.txt",
            "salmonrun/sleuth_completed.txt",
            expand("{name}.RnaSeqMetrics.txt",name=samples),
            "postTrimQC","sampletable.txt",
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
            "RawCountFileStar.txt",expand("{name}.rsem.genes.results",name=samples)

elif config['project']['DEG'] == "no" and config['project']['TRIM'] == "yes":
  rule all:
     params: batch='--time=168:00:00'
     input: "STAR_QC","Reports/multiqc_report.html",
            expand("{name}.RnaSeqMetrics.txt",name=samples),
            "postTrimQC",
            "RawCountFile_genes_filtered.txt",
            "RawCountFile_junctions_filtered.txt",
            "RawCountFile_genejunctions_filtered.txt",
            expand("{name}.star.count.overlap.txt",name=samples),
            "RawCountFileOverlap.txt","RawCountFileStar.txt",
            expand("{name}.rsem.genes.results",name=samples)

elif config['project']['DEG'] == "yes" and config['project']['TRIM'] == "no":
  rule all:
     input: "STAR_QC","Reports/multiqc_report.html",
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
            expand("{name}.rsem.genes.results",name=samples)
            
     params: batch='--time=168:00:00'
     #input: "files_to_rnaseqc.txt","STAR_QC","RawCountFile_filtered.txt","sampletable.txt","deseq2_pca.png","edgeR_prcomp.png","Limma_MDS.png"
else:
  rule all:
     params: batch='--time=168:00:00'
     input: "STAR_QC","Reports/multiqc_report.html",
            expand("{name}.RnaSeqMetrics.txt",name=samples),
            "RawCountFile_genes_filtered.txt",
            "RawCountFile_genejunctions_filtered.txt",
            expand("{name}.star.count.overlap.txt",name=samples),
            "RawCountFileOverlap.txt","RawCountFileStar.txt",
            expand("{name}.rsem.genes.results",name=samples)


if config['project']['TRIM'] == "yes":
   rule trimmomatic_pe:
      input: file1= config['project']['workpath']+"/{name}.R1."+config['project']['filetype'],file2=config['project']['workpath']+"/{name}.R2."+config['project']['filetype'] 
      output: out11="trim/{name}_R1_001_trim_paired.fastq.gz",out12="trim/{name}_R1_001_trim_unpaired.fastq.gz",out21="trim/{name}_R2_001_trim_paired.fastq.gz",out22="trim/{name}_R2_001_trim_unpaired.fastq.gz"
      params: rname='pl:trimmomatic_pe',batch='--cpus-per-task=32 --mem=110g --time=48:00:00',trimmomaticver=config['bin'][pfamily]['TRIMMOMATICVER'],fastawithadaptersetc=config['references'][pfamily]['FASTAWITHADAPTERSETC'],seedmismatches=config['bin'][pfamily]['SEEDMISMATCHES'],palindromeclipthreshold=config['bin'][pfamily]['PALINDROMECLIPTHRESHOLD'],simpleclipthreshold=config['bin'][pfamily]['SIMPLECLIPTHRESHOLD'],leadingquality=config['bin'][pfamily]['LEADINGQUALITY'],trailingquality=config['bin'][pfamily]['TRAILINGQUALITY'],windowsize=config['bin'][pfamily]['WINDOWSIZE'],windowquality=config['bin'][pfamily]['WINDOWQUALITY'],targetlength=config['bin'][pfamily]['TARGETLENGTH'],strictness=config['bin'][pfamily]['STRICTNESS'],minlen=config['bin'][pfamily]['MINLEN'],headcroplength=config['bin'][pfamily]['HEADCROPLENGTH']
      threads:32
      shell:"module load {params.trimmomaticver}; java -classpath $TRIMMOJAR   org.usadellab.trimmomatic.TrimmomaticPE -threads {threads} {input.file1} {input.file2} {output.out11} {output.out12} {output.out21} {output.out22} ILLUMINACLIP:{params.fastawithadaptersetc}:{params.seedmismatches}:{params.palindromeclipthreshold}:{params.simpleclipthreshold}  LEADING:{params.leadingquality} TRAILING:{params.trailingquality} SLIDINGWINDOW:{params.windowsize}:{params.windowquality} MAXINFO:{params.targetlength}:{params.strictness} MINLEN:{params.minlen} HEADCROP:{params.headcroplength}"

# ILLUMINACLIP:{params.fastawithadaptersetc}:{params.seedmismatches}:{params.palindromeclipthreshold}:{params.simpleclipthreshold}  LEADING:{params.leadingquality} TRAILING:{params.trailingquality} SLIDINGWINDOW:{params.windowsize}:{params.windowquality} MAXINFO:{params.targetlength}:{params.strictness} MINLEN:{params.minlen} HEADCROP:{params.headcroplength}"


   rule fastqc:  
      input: expand("trim/{name}_R1_001_trim_paired.fastq.gz", name=samples), expand("trim/{name}_R2_001_trim_paired.fastq.gz", name=samples)  
      output: "postTrimQC"
      priority: 2
      params: rname='pl:fastqc',batch='--cpus-per-task=32 --mem=110g --time=48:00:00',fastqcver=config['bin'][pfamily]['FASTQCVER']
      threads: 32
      shell: "mkdir -p {output};module load {params.fastqcver}; fastqc {input} -t {threads} -o {output}"

#   rule check:
#      input: "postTrimQC"
#      output: out1="FastqcSummary.xlsx", out2="fastqc_status.txt"
#      params: batch='-l nodes=1:gpfs',chkqcver=config['bin'][pfamily]['CHKQCVER']
#      shell: """
#             {params.chkqcver} --dir {input} --outfile {output.out1};
#             if [ "$?" == 1 ]; then 
#                 echo "FAILURE" > {output.out2}
#             else 
#                 echo "OK" > {output.out2} 
#             fi
#             """
#   rule decide:
#      input: file1="fastqc_status.txt"
#      output: out1="fastqc_status_checked.txt"
#      params: batch='--time=4:00:00',resume=config['bin'][pfamily]['RESUME']
#      shell: """
#             status=$(<{input.file1})
#             if [ $status == "OK" -o {params.resume} == "yes" ]; then
#                  if [ {params.resume} == "yes" -a $status == "FAILURE"  ]; then
#                      echo "workflow resumed after QC check" > {output.out1}
#                  else
#                      echo "status OK" > {output.out1}
#                  fi 
#             else
#                  exit 1
#             fi
#            """

   rule star1p:
      input: file1= "trim/{name}_R1_001_trim_paired.fastq.gz",file2="trim/{name}_R2_001_trim_paired.fastq.gz"#,file3="fastqc_status_checked.txt"
      output: out1= "{name}.SJ.out.tab"#,out2= "{name}.SJ.out.tab.Pass1.sjdb"
      params: rname='pl:star1p',prefix="{name}",batch='--cpus-per-task=32 --mem=110g --time=48:00:00',starver=config['bin'][pfamily]['STARVER'],stardir=config['references'][pfamily]['STARDIR']+config['project']["SJDBOVERHANG"],filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],filtertype=config['bin'][pfamily]['FILTERTYPE'],filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],adapter1=config['bin'][pfamily]['ADAPTER1'],adapter2=config['bin'][pfamily]['ADAPTER2'],
      threads: 32
      shell: """
             module load {params.starver}
             STAR --genomeDir {params.stardir} --outFilterIntronMotifs {params.filterintronmotifs} --outSAMstrandField {params.samstrandfield}  --outFilterType {params.filtertype} --outFilterMultimapNmax {params.filtermultimapnmax} --alignSJoverhangMin {params.alignsjoverhangmin} --alignSJDBoverhangMin {params.alignsjdboverhangmin}  --outFilterMismatchNmax {params.filtermismatchnmax} --outFilterMismatchNoverLmax {params.filtermismatchnoverlmax}  --alignIntronMin {params.alignintronmin} --alignIntronMax {params.alignintronmax} --alignMatesGapMax {params.alignmatesgapmax} clip3pAdapterSeq {params.adapter1} {params.adapter2} --readFilesIn {input.file1} {input.file2} --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {params.prefix}. --outSAMtype BAM Unsorted
             """  
#             awk 'BEGIN {{OFS=\"\\t\"; strChar[0]=\".\"; strChar[1]=\"+\"; strChar[2]=\"-\";}} {{if($5>0){{print $1,$2,$3,strChar[$4]}}}}' {output.out1}  > {output.out2}



   rule star2p:
      input: file1= "trim/{name}_R1_001_trim_paired.fastq.gz",file2="trim/{name}_R2_001_trim_paired.fastq.gz",tab=expand("{name}.SJ.out.tab",name=samples)#, dir="STARINDEX"
      output: out1="{name}.p2.Aligned.sortedByCoord.out.bam",out2="{name}.p2.ReadsPerGene.out.tab",out3="{name}.p2.Aligned.toTranscriptome.out.bam",out4="{name}.p2.SJ.out.tab" #"{name}.p2.Aligned.out.sam"
      params: rname='pl:star2p',prefix="{name}.p2",batch='--cpus-per-task=32 --mem=110g --time=48:00:00',starver=config['bin'][pfamily]['STARVER'],filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],filtertype=config['bin'][pfamily]['FILTERTYPE'],filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],adapter1=config['bin'][pfamily]['ADAPTER1'],adapter2=config['bin'][pfamily]['ADAPTER2'],outsamunmapped=config['bin'][pfamily]['OUTSAMUNMAPPED'],wigtype=config['bin'][pfamily]['WIGTYPE'],wigstrand=config['bin'][pfamily]['WIGSTRAND'], gtffile=config['references'][pfamily]['GTFFILE'], nbjuncs=config['bin'][pfamily]['NBJUNCS'],stardir=config['references'][pfamily]['STARDIR']+config['project']["SJDBOVERHANG"]
      threads:32
      shell:"module load {params.starver}; STAR --genomeDir {params.stardir} --outFilterIntronMotifs {params.filterintronmotifs} --outSAMstrandField {params.samstrandfield}  --outFilterType {params.filtertype} --outFilterMultimapNmax {params.filtermultimapnmax} --alignSJoverhangMin {params.alignsjoverhangmin} --alignSJDBoverhangMin {params.alignsjdboverhangmin}  --outFilterMismatchNmax {params.filtermismatchnmax} --outFilterMismatchNoverLmax {params.filtermismatchnoverlmax}  --alignIntronMin {params.alignintronmin} --alignIntronMax {params.alignintronmax} --alignMatesGapMax {params.alignmatesgapmax}  --clip3pAdapterSeq {params.adapter1} {params.adapter2} --readFilesIn {input.file1} {input.file2} --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {params.prefix}. --outSAMunmapped {params.outsamunmapped} --outWigType {params.wigtype} --outWigStrand {params.wigstrand} --sjdbFileChrStartEnd {input.tab} --sjdbGTFfile {params.gtffile} --limitSjdbInsertNsj {params.nbjuncs} --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate"


else:
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



#rule newindex:
#  input: expand("{name}.SJ.out.tab.Pass1.sjdb", name=samples)
#  output: dir="STARINDEX",file="all.SJ.out.tab.Pass1.sjdb"
#  params: batch='--cpus-per-task=32 --mem=110g --time=48:00:00',starver=config['bin'][pfamily]['STARVER'],genomefile=config['references'][pfamily]['GENOMEFILE'],gtffile=config['references'][pfamily]['GTFFILE'],sjdboverhang=config['bin'][pfamily]['SJDBOVERHANG']
#  threads: 32
#  shell:"mkdir {output.dir}; cat {input} > {output.file}; module load {params.starver}; STAR --runMode genomeGenerate --genomeDir {output.dir} --genomeFastaFiles {params.genomefile} --sjdbGTFfile {params.gtffile} --sjdbFileChrStartEnd {output.file} --sjdbOverhang {params.sjdboverhang} --runThreadN {threads}" 

rule rsem:
  input: file1= "{name}.p2.Aligned.toTranscriptome.out.bam"
  output: out1="{name}.rsem.genes.results",out2="{name}.rsem.isoforms.results"
  params: rname='pl:rsem',prefix="{name}.rsem",batch='--cpus-per-task=16 --mem=32g --time=24:00:00',rsemref=config['references'][pfamily]['RSEMREF'],rsem=config['bin'][pfamily]['RSEM']
  shell: "{params.rsem}/rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  --bam --paired-end -p 16  {input.file1} {params.rsemref} {params.prefix}"

rule picard:
  input: file1= "{name}.p2.Aligned.sortedByCoord.out.bam"
  output: outstar1="{name}.star_rg_added.sorted.bam", outstar2="{name}.star_rg_added.sorted.dmark.bam",outstar3="{name}.star.duplic" 
  params: rname='pl:picard',batch='--mem=24g --time=10:00:00 --gres=lscratch:800',picardver=config['bin'][pfamily]['PICARDVER']#,picardjarpath=config['bin'][pfamily]['PICARDJARPATH']
  shell: "module load {params.picardver}; java -Xmx10g  -jar $PICARDJARPATH/AddOrReplaceReadGroups.jar INPUT={input.file1} OUTPUT={output.outstar1} TMP_DIR=/lscratch/$SLURM_JOBID RGID=id RGLB=library RGPL=illumina RGPU=machine RGSM=sample; java -Xmx10g -jar $PICARDJARPATH/MarkDuplicates.jar INPUT={output.outstar1} OUTPUT={output.outstar2} TMP_DIR=/lscratch/$SLURM_JOBID CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE={output.outstar3}"

rule stats:
  input: file1= "{name}.star_rg_added.sorted.dmark.bam"
  output: outstar1="{name}.RnaSeqMetrics.txt", outstar2="{name}.flagstat.concord.txt"
  params: rname='pl:stats',batch='--mem=24g --time=10:00:00 --gres=lscratch:800',picardver=config['bin'][pfamily]['PICARDVER'],refflat=config['references'][pfamily]['REFFLAT'],rrnalist=config['references'][pfamily]['RRNALIST'],picardstrand=config['bin'][pfamily]['PICARDSTRAND']
  shell: "module load {params.picardver}; java -Xmx10g -jar $PICARDJARPATH/CollectRnaSeqMetrics.jar REF_FLAT={params.refflat} INPUT={input.file1} OUTPUT={output.outstar1} RIBOSOMAL_INTERVALS={params.rrnalist}  STRAND_SPECIFICITY={params.picardstrand} TMP_DIR=/lscratch/$SLURM_JOBID  VALIDATION_STRINGENCY=SILENT; module load samtools; samtools flagstat {input.file1} > {output.outstar2}; samtools  view -f 0x2 {input.file1} | wc -l >>{output.outstar2}; samtools view {input.file1} | grep -w -c NH:i:1  >>{output.outstar2} "

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

rule subread:
   input:  "{name}.star_rg_added.sorted.dmark.bam"
   output: out="{name}.star.count.info.txt", res="{name}.star.count.txt"
   params: rname='pl:subread',batch='--time=4:00:00',subreadver=config['bin'][pfamily]['SUBREADVER'],stranded=config['bin'][pfamily]['STRANDED'],gtffile=config['references'][pfamily]['GTFFILE']
   shell: "module load {params.subreadver}; featureCounts -T 16 -s {params.stranded} -p -t exon -R -g gene_id -a {params.gtffile} -o {output.out}  {input}; sed '1d' {output.out} | cut -f1,7 > {output.res}"

rule subreadoverlap:
   input:  "{name}.star_rg_added.sorted.dmark.bam"
   output: out="{name}.star.count.info.overlap.txt", res="{name}.star.count.overlap.txt"
   params: rname='pl:subreadoverlap',batch='--cpus-per-task=16 --mem=24g --time=48:00:00',subreadver=config['bin'][pfamily]['SUBREADVER'],stranded=config['bin'][pfamily]['STRANDED'],gtffile=config['references'][pfamily]['GTFFILE']
   shell: "module load {params.subreadver}; featureCounts -T 16 -s {params.stranded} -p -t exon -R -O -g gene_id -a {params.gtffile} -o {output.out}  {input}; sed '1d' {output.out} | cut -f1,7 > {output.res}"


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
    input: "STAR_QC/index.html","STAR_QC/report.html"
    output: "Reports/multiqc_report.html"
    params: rname="pl:multiqc",pythonpath=config['bin'][pfamily]['PYTHONPATH'],multiqc=config['bin'][pfamily]['MULTIQC']
    threads: 1
    shell:  """
            module load multiqc
            cd Reports && multiqc -f -e featureCounts -e picard ../

            """
   
rule samplecondition:
   input: files=expand("{name}.star.count.txt", name=samples)
   output: out1= "sampletable.txt"
   params: rname='pl:samplecondition',batch='--mem=4g --time=10:00:00', groups=config['project']['contrasts']['rgroups'], labels=config['project']['contrasts']['rlabels']
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
  shell: "mkdir -p DEG{params.dtype}; module load R; Rscript Scripts/deseq2.R '{params.dir}/DEG{params.dtype}/' '../{input.file1}' '../{input.file2}' '{params.annotate}' '{params.contrasts}'"

rule edgeR:
  input: file1="sampletable.txt", file2="RawCountFile{dtype}_filtered.txt"
  output: "DEG{dtype}/edgeR_prcomp.png"
  params: rname='pl:edgeR',batch='--mem=24g --time=10:00:00', dir=config['project']['workpath'],annotate=config['references'][pfamily]['ANNOTATE'],contrasts=" ".join(config['project']['contrasts']['rcontrasts']), dtype="{dtype}"
  shell: "mkdir -p DEG{params.dtype}; module load R; Rscript Scripts/edgeR.R '{params.dir}/DEG{params.dtype}/' '../{input.file1}' '../{input.file2}' '{params.annotate}' '{params.contrasts}'"

rule limmavoom:
  input: file1="sampletable.txt", file2="RawCountFile{dtype}_filtered.txt"
  output: "DEG{dtype}/Limma_MDS.png"
  params: rname='pl:limmavoom',batch='--mem=24g --time=10:00:00',dir=config['project']['workpath'],contrasts=" ".join(config['project']['contrasts']['rcontrasts']), dtype="{dtype}"
  shell: "mkdir -p DEG{params.dtype}; module load R; Rscript Scripts/limmavoom.R '{params.dir}/DEG{params.dtype}/' '../{input.file1}' '../{input.file2}' '{params.contrasts}'"

rule salmon:
  input: bam="{name}.p2.Aligned.toTranscriptome.out.bam"
  output: "salmonrun/{name}/quant.sf"
  params: sname="{name}",rname='pl:salmon',batch='--mem=128g --cpus-perptask=8 --time=10:00:00',dir=config['project']['workpath'],rsemref=config['references'][pfamily]['SALMONREF'],libtype={0:'U',1:'SF',2:'SR'}.get(config['bin'][pfamily]['STRANDED'])
  shell: "mkdir -p {params.dir}/salmonrun; module load salmon/0.6.0; salmon quant -t {params.rsemref} -l I{params.libtype} -a {input.bam} -o {params.dir}/salmonrun/{params.sname} --numBootstraps 30;"

rule sleuth:
  input: samtab = "sampletable.txt", bam=expand("salmonrun/{name}/quant.sf", name=samples)
  output: "salmonrun/sleuth_completed.txt"
  params: rname='pl:sleuth',batch='--mem=128g --cpus-per-task=8 --time=10:00:00',dir=config['project']['workpath'],pipeRlib=config['bin'][pfamily]['PIPERLIB'],contrasts=" ".join(config['project']['contrasts']['rcontrasts'])
  shell: "module load R; Rscript Scripts/sleuth.R '{params.dir}' '{params.pipeRlib}' '{input.samtab}' '{params.contrasts}'"

rule EBSeq:
  input: samtab = "sampletable.txt", isoforms=expand("{name}.rsem.isoforms.results", name=samples)
  output: "ebseq_completed.txt"
  params: rname='pl:EBSeq',batch='--mem=128g --cpus-per-task=8 --time=10:00:00',dir=config['project']['workpath'],contrasts=" ".join(config['project']['contrasts']['rcontrasts']),rsemref=config['references'][pfamily]['RSEMREF'],rsem=config['bin'][pfamily]['RSEM']
  shell: "module load R; Rscript Scripts/ebseq.R '{params.dir}' '{input.isoforms}' '{input.samtab}' '{params.contrasts}' '{params.rsemref}' '{params.rsem}'"
