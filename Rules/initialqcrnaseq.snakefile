from snakemake.utils import R
from os.path import join
from os import listdir

configfile: "run.json"
samples=config['project']['groups']['rsamps']

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

se=""
pe=""
workpath = config['project']['workpath']

if config['project']['nends'] == 2 :
	pe="yes"
elif config['project']['nends'] == 1 :
	se="yes"


trim_dir='trim'
star_dir="STAR_files"
bams_dir="bams"
log_dir="logfiles"
rseqc_dir="RSeQC"
kraken_dir="kraken"
preseq_dir="preseq"
degall_dir="DEG_ALL"

dtypes=["RSEM_genes","Subread_junctions","Subread_genejunctions","Subread_genes"]

for d in [trim_dir,kraken_dir,bams_dir,star_dir,log_dir,rseqc_dir,preseq_dir,degall_dir]:
	if not os.path.exists(join(workpath,d)):
		os.mkdir(join(workpath,d))


if pe=="yes":

   rule all:
      params:
        batch='--time=168:00:00',
        # rname='pl:all',
      input:
        # FastQC (before and after trimming)
        join(workpath,"rawQC"),
        join(workpath,"QC"),

        # FastScreen (Using two sets of reference databases)
        expand(join(workpath,"FQscreen","{name}.R1.trim_screen.txt"),name=samples),
        expand(join(workpath,"FQscreen","{name}.R1.trim_screen.png"),name=samples),
        expand(join(workpath,"FQscreen","{name}.R2.trim_screen.txt"),name=samples),
        expand(join(workpath,"FQscreen","{name}.R2.trim_screen.png"),name=samples),
        expand(join(workpath,"FQscreen2","{name}.R1.trim_screen.txt"),name=samples),
        expand(join(workpath,"FQscreen2","{name}.R1.trim_screen.png"),name=samples),
        expand(join(workpath,"FQscreen2","{name}.R2.trim_screen.txt"),name=samples),
        expand(join(workpath,"FQscreen2","{name}.R2.trim_screen.png"),name=samples),

        # Kraken + Krona
        expand(join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.taxa.txt"),name=samples),
        expand(join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.krona.html"),name=samples),

        # STAR
        expand(join(workpath,bams_dir,"{name}.p2.Aligned.toTranscriptome.out.bam"),name=samples),

        # Picard
        expand(join(workpath,log_dir,"{name}.RnaSeqMetrics.txt"),name=samples),
        expand(join(workpath,log_dir,"{name}.star.duplic"),name=samples),

        # Preseq
        expand(join(workpath,preseq_dir,"{name}.ccurve"),name=samples),

        # Infer Library Strandedness
        join(workpath,log_dir,"strandness.txt"),

        # QualiMap (bamQC and counts)
        expand(join(workpath,"QualiMap","{name}","qualimapReport.html"),name=samples),
        # join(workpath,"QualiMap","GlobalReport.html"),

        # RSEM merge and counts 
        join(workpath,degall_dir,"RSEM.genes.FPKM.all_samples.txt"),
        join(workpath,degall_dir,"RSEM.isoforms.FPKM.all_samples.txt"),
        join(workpath,degall_dir,"RawCountFile_RSEM_genes_filtered.txt"),

        # Subread 
        expand(join(workpath,star_dir,"{name}.star.count.info.txt"),name=samples),
        expand(join(workpath,star_dir,"{name}.star.count.txt"),name=samples),

        # Subread featureCounts
        join(workpath,star_dir,"sampletable.txt"),
        join(workpath,degall_dir,"RawCountFile_Subread_genes_filtered.txt"),
        join(workpath,degall_dir,"RawCountFile_Subread_junctions_filtered.txt"),
        join(workpath,degall_dir,"RawCountFile_Subread_genejunctions_filtered.txt"),

        # Subread Overlap 
        join(workpath,star_dir,"RawCountFileOverlap.txt"),
        join(workpath,star_dir,"RawCountFileStar.txt"),    	
        expand(join(workpath,star_dir,"{name}.star.count.info.overlap.txt"),name=samples),
        expand(join(workpath,star_dir,"{name}.star.count.overlap.txt"),name=samples),

        # PCA Reports 
        expand(join(workpath,degall_dir,"PcaReport_{dtype}.html"),dtype=dtypes),
        
        # MultiQC
        join(workpath,"Reports/multiqc_report.html"),

   rule rawfastqc:
      input: 
        expand("{name}.R1.fastq.gz", name=samples), 
        expand("{name}.R2.fastq.gz", name=samples)
      output: 
        join(workpath,"rawQC")
      priority: 2
      params: 
        rname='pl:rawfastqc',
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        fastqcver=config['bin'][pfamily]['tool_versions']['FASTQCVER'],
      threads: 32
      shell: """
mkdir -p {output};
module load {params.fastqcver};
fastqc {input} -t {threads} -o {output}; 
        """

   rule trim_pe:
      input: 
        file1=join(workpath,"{name}.R1."+config['project']['filetype']),
        file2=join(workpath,"{name}.R2."+config['project']['filetype']),        
      output: 
        out1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        out2=join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
      params: 
        rname='pl:trim_pe',
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        cutadaptver=config['bin'][pfamily]['tool_versions']['CUTADAPTVER'],
        fastawithadaptersetd=config['bin'][pfamily]['tool_parameters']['FASTAWITHADAPTERSETD'],
        leadingquality=config['bin'][pfamily]['tool_parameters']['LEADINGQUALITY'],
        trailingquality=config['bin'][pfamily]['tool_parameters']['TRAILINGQUALITY'],
        minlen=config['bin'][pfamily]['tool_parameters']['MINLEN'],
      threads:32
      shell: """
module load {params.cutadaptver};
sample=`echo {input.file1}|awk -F "/" '{{print $NF}}'|awk -F ".R1.fastq" '{{print $1}}'`
cutadapt --pair-filter=any --nextseq-trim=2 --trim-n -n 5 -O 5 -q {params.leadingquality},{params.trailingquality} -m {params.minlen}:{params.minlen} -b file:{params.fastawithadaptersetd} -B file:{params.fastawithadaptersetd} -j {threads} -o {output.out1} -p {output.out2} {input.file1} {input.file2}
        """

   rule fastqc:
      input: 
        expand(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"), name=samples), 
        expand(join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"), name=samples)
      output: 
        join(workpath,"QC")
      priority: 2
      params: 
        rname='pl:fastqc',
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        fastqcver=config['bin'][pfamily]['tool_versions']['FASTQCVER'],
      threads: 32
      shell: """
mkdir -p {output};
module load {params.fastqcver};
fastqc {input} -t {threads} -o {output};
module load python/3.5;
python Scripts/get_read_length.py {output} > {output}/readlength.txt  2> {output}/readlength.err
        """

   rule fastq_screen:
      input: 
        file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"), 
        file2=join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
      output: 
        out1=join(workpath,"FQscreen","{name}.R1.trim_screen.txt"), 
        out2=join(workpath,"FQscreen","{name}.R1.trim_screen.png"), 
        out3=join(workpath,"FQscreen","{name}.R2.trim_screen.txt"), 
        out4=join(workpath,"FQscreen","{name}.R2.trim_screen.png"),
        out5=join(workpath,"FQscreen2","{name}.R1.trim_screen.txt"), 
        out6=join(workpath,"FQscreen2","{name}.R1.trim_screen.png"), 
        out7=join(workpath,"FQscreen2","{name}.R2.trim_screen.txt"), 
        out8=join(workpath,"FQscreen2","{name}.R2.trim_screen.png")
      params: 
        rname='pl:fqscreen',
        batch='--cpus-per-task=24 --mem=64g --time=10:00:00',
        fastq_screen=config['bin'][pfamily]['tool_versions']['FASTQ_SCREEN'],
        outdir = join(workpath,"FQscreen"),
        outdir2 = join(workpath,"FQscreen2"),
        fastq_screen_config=config['bin'][pfamily]['tool_parameters']['FASTQ_SCREEN_CONFIG'],
        fastq_screen_config2=config['bin'][pfamily]['tool_parameters']['FASTQ_SCREEN_CONFIG2'],
        perlver=config['bin'][pfamily]['tool_versions']['PERLVER'],
        bowtie2ver=config['bin'][pfamily]['tool_versions']['BOWTIE2VER'],
      threads: 24
      shell: """
module load {params.bowtie2ver};
module load {params.perlver};
{params.fastq_screen} --conf {params.fastq_screen_config} --outdir {params.outdir} --threads {threads} --subset 1000000 --aligner bowtie2 --force {input.file1} {input.file2}
{params.fastq_screen} --conf {params.fastq_screen_config2} --outdir {params.outdir2} --threads {threads} --subset 1000000 --aligner bowtie2 --force {input.file1} {input.file2}
        """

   rule kraken_pe:
      input: 
        fq1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        fq2=join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
      output: 
        krakentaxa = join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.taxa.txt"),
        kronahtml = join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.krona.html"),
      params: 
        rname='pl:kraken',
        prefix = "{name}",
        outdir=join(workpath,kraken_dir),
        bacdb=config['bin'][pfamily]['tool_parameters']['KRAKENBACDB'],
        krakenver=config['bin'][pfamily]['tool_versions']['KRAKENVER'],
        kronatoolsver=config['bin'][pfamily]['tool_versions']['KRONATOOLSVER'],
      threads: 24
      shell: """
module load {params.krakenver};
module load {params.kronatoolsver};
if [ ! -d {params.outdir} ];then mkdir {params.outdir};fi
cd /lscratch/$SLURM_JOBID;
cp -rv {params.bacdb} /lscratch/$SLURM_JOBID/;
kraken --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` --fastq-input --gzip-compressed --threads {threads} --output /lscratch/$SLURM_JOBID/{params.prefix}.krakenout --preload --paired {input.fq1} {input.fq2}
kraken-translate --mpa-format --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` /lscratch/$SLURM_JOBID/{params.prefix}.krakenout |cut -f2|sort|uniq -c|sort -k1,1nr > /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa
cut -f2,3 /lscratch/$SLURM_JOBID/{params.prefix}.krakenout | ktImportTaxonomy - -o /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml
mv /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa {output.krakentaxa}
mv /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml {output.kronahtml}
"""


   rule star1p:
      input: 
        file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        file2=join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
        qcdir=join(workpath,"QC"),
        # file3="fastqc_status_checked.txt"
      output: 
        out1= join(workpath,star_dir,"{name}.SJ.out.tab"), 
        out3= temp(join(workpath,star_dir,"{name}.Aligned.out.bam")),
        # out2= "{name}.SJ.out.tab.Pass1.sjdb"
      params: 
        rname='pl:star1p',prefix="{name}",
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        starver=config['bin'][pfamily]['tool_versions']['STARVER'],
        stardir=config['references'][pfamily]['STARDIR'],
        filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],
        samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],
        filtertype=config['bin'][pfamily]['FILTERTYPE'],
        filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],
        alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],
        alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],
        filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],
        filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],
        alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],
        alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],
        alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],
        adapter1=config['bin'][pfamily]['ADAPTER1'],
        adapter2=config['bin'][pfamily]['ADAPTER2'],
      threads: 32
      run:
        import glob,json
        rl=int(open(join(input.qcdir,"readlength.txt")).readlines()[0].strip())-1
        dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.stardir+'*/',recursive=False))))
        bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
        cmd="cd {star_dir}; module load "+params.starver+"; STAR --genomeDir "+params.stardir+str(bestdbrl)+" --outFilterIntronMotifs "+params.filterintronmotifs+" --outSAMstrandField "+params.samstrandfield+"  --outFilterType "+params.filtertype+" --outFilterMultimapNmax "+str(params.filtermultimapnmax)+" --alignSJoverhangMin "+str(params.alignsjoverhangmin)+" --alignSJDBoverhangMin "+str(params.alignsjdboverhangmin)+"  --outFilterMismatchNmax "+str(params.filtermismatchnmax)+" --outFilterMismatchNoverLmax "+str(params.filtermismatchnoverlmax)+"  --alignIntronMin "+str(params.alignintronmin)+" --alignIntronMax "+str(params.alignintronmax)+" --alignMatesGapMax "+str(params.alignmatesgapmax)+" --clip3pAdapterSeq "+params.adapter1+" "+params.adapter2+" --readFilesIn "+input.file1+" "+input.file2+" --readFilesCommand zcat --runThreadN "+str(threads)+" --outFileNamePrefix "+params.prefix+". --outSAMtype BAM Unsorted"
        shell(cmd)
        A=open(join(workpath,"run.json"),'r')
        a=eval(A.read())
        A.close()
        config=dict(a.items())
        config['project']['SJDBOVERHANG']=str(bestdbrl)
        config['project']["STARDIR"]= config['references'][pfamily]['STARDIR']+str(bestdbrl)
        config['project']['READLENGTH']=str(rl+1)
        with open(join(workpath,'run.json'),'w') as F:
          json.dump(config, F, sort_keys = True, indent = 4,ensure_ascii=False)
        F.close()

   rule sjdb:
      input: 
        files=expand(join(workpath,star_dir,"{name}.SJ.out.tab"), name=samples)
      output: 
        out1=join(workpath,star_dir,"uniq.filtered.SJ.out.tab")
      params: 
        rname='pl:sjdb'
      shell: """
cat {input.files} |sort|uniq|awk -F \"\\t\" '{{if ($5>0 && $6==1) {{print}}}}'|cut -f1-4|sort|uniq|grep \"^chr\"|grep -v \"^chrM\" > {output.out1}
        """

   rule star2p:
      input: 
        file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        file2=join(workpath,trim_dir,"{name}.R2.trim.fastq.gz"),
        tab=join(workpath,star_dir,"uniq.filtered.SJ.out.tab"),
        qcdir=join(workpath,"QC"),
        # dir="STARINDEX"
      output: 
        out1=temp(join(workpath,star_dir,"{name}.p2.Aligned.sortedByCoord.out.bam")),
        out2=join(workpath,star_dir,"{name}.p2.ReadsPerGene.out.tab"),
        out3=join(workpath,bams_dir,"{name}.p2.Aligned.toTranscriptome.out.bam"),
        out4=join(workpath,star_dir,"{name}.p2.SJ.out.tab"),
        out5=join(workpath,log_dir,"{name}.p2.Log.final.out"),
      params: 
        rname='pl:star2p',
        prefix="{name}.p2",
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        starver=config['bin'][pfamily]['tool_versions']['STARVER'],
        filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],
        samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],
        filtertype=config['bin'][pfamily]['FILTERTYPE'],
        filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],
        alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],
        alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],
        filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],
        filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],
        alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],
        alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],
        alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],
        adapter1=config['bin'][pfamily]['ADAPTER1'],
        adapter2=config['bin'][pfamily]['ADAPTER2'],
        outsamunmapped=config['bin'][pfamily]['OUTSAMUNMAPPED'],
        wigtype=config['bin'][pfamily]['WIGTYPE'],
        wigstrand=config['bin'][pfamily]['WIGSTRAND'],
        gtffile=config['references'][pfamily]['GTFFILE'],
        nbjuncs=config['bin'][pfamily]['NBJUNCS'],
        stardir=config['references'][pfamily]['STARDIR'],
      threads:32
      run:
        import glob,os
        rl=int(open(join(input.qcdir,"readlength.txt")).readlines()[0].strip())-1
        dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.stardir+'*/',recursive=False))))
        bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
        cmd="cd {star_dir}; module load "+params.starver+"; STAR --genomeDir "+params.stardir+str(bestdbrl)+" --outFilterIntronMotifs "+params.filterintronmotifs+" --outSAMstrandField "+params.samstrandfield+"  --outFilterType "+params.filtertype+" --outFilterMultimapNmax "+str(params.filtermultimapnmax)+" --alignSJoverhangMin "+str(params.alignsjoverhangmin)+" --alignSJDBoverhangMin "+str(params.alignsjdboverhangmin)+"  --outFilterMismatchNmax "+str(params.filtermismatchnmax)+" --outFilterMismatchNoverLmax "+str(params.filtermismatchnoverlmax)+"  --alignIntronMin "+str(params.alignintronmin)+" --alignIntronMax "+str(params.alignintronmax)+" --alignMatesGapMax "+str(params.alignmatesgapmax)+"  --clip3pAdapterSeq "+params.adapter1+" "+params.adapter2+" --readFilesIn "+input.file1+" "+input.file2+" --readFilesCommand zcat --runThreadN "+str(threads)+" --outFileNamePrefix "+params.prefix+". --outSAMunmapped "+params.outsamunmapped+" --outWigType "+params.wigtype+" --outWigStrand "+params.wigstrand+" --sjdbFileChrStartEnd "+str(input.tab)+" --sjdbGTFfile "+params.gtffile+" --limitSjdbInsertNsj "+str(params.nbjuncs)+" --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate;"
        shell(cmd)
        cmd="sleep 120;cd {workpath};mv {workpath}/{star_dir}/{params.prefix}.Aligned.toTranscriptome.out.bam {workpath}/{bams_dir}; mv {workpath}/{star_dir}/{params.prefix}.Log.final.out {workpath}/{log_dir}"
        shell(cmd)

   rule qualibam:
      input:
        bamfile=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
      output:
        report=join(workpath,"QualiMap","{name}","qualimapReport.html"),
      params:
        rname='pl:qualibam',
        outdir=join(workpath,"QualiMap","{name}"),
        gtfFile=config['references'][pfamily]['GTFFILE'],
      shell: """
module load qualimap/2.2.1
unset DISPLAY;qualimap bamqc -bam {input.bamfile} --feature-file {params.gtfFile} -outdir {params.outdir} -nt $SLURM_CPUS_PER_TASK --java-mem-size=11G
        """


if se=="yes":

   rule all:
      params: batch='--time=168:00:00'
      input:
        # FastQC (before and after trimming)
        join(workpath,"rawQC"),
        join(workpath,"QC"),

        # FastQ Screen
        expand(join(workpath,"FQscreen","{name}.R1.trim_screen.txt"),name=samples),
        expand(join(workpath,"FQscreen","{name}.R1.trim_screen.png"),name=samples),
        expand(join(workpath,"FQscreen2","{name}.R1.trim_screen.txt"),name=samples),
        expand(join(workpath,"FQscreen2","{name}.R1.trim_screen.png"),name=samples),

        # Kraken + Krona
        expand(join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.taxa.txt"),name=samples),
        expand(join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.krona.html"),name=samples),

        # Infer Library Strandedness
        join(workpath,log_dir,"strandness.txt"),

        # Picard 
        expand(join(workpath,log_dir,"{name}.RnaSeqMetrics.txt"),name=samples),

        # QualiMap
        # join(workpath,"QualiMap","GlobalReport.html"),

        # RSEM
        join(workpath,degall_dir,"RSEM.genes.FPKM.all_samples.txt"),
        join(workpath,degall_dir,"RSEM.isoforms.FPKM.all_samples.txt"),
        join(workpath,degall_dir,"RawCountFile_RSEM_genes_filtered.txt"),

        # Subread
        expand(join(workpath,star_dir,"{name}.star.count.info.txt"),name=samples),
        expand(join(workpath,star_dir,"{name}.star.count.txt"),name=samples),

        # Subread featureCounts
        join(workpath,star_dir,"sampletable.txt"),
        join(workpath,degall_dir,"RawCountFile_Subread_genes_filtered.txt"),
        join(workpath,degall_dir,"RawCountFile_Subread_junctions_filtered.txt"),
        join(workpath,degall_dir,"RawCountFile_Subread_genejunctions_filtered.txt"),

        # Subread Overlaps
        join(workpath,star_dir,"RawCountFileOverlap.txt"),
        join(workpath,star_dir,"RawCountFileStar.txt"),
        expand(join(workpath,star_dir,"{name}.star.count.info.overlap.txt"),name=samples),
        expand(join(workpath,star_dir,"{name}.star.count.overlap.txt"),name=samples),

        # PCA Report
        expand(join(workpath,degall_dir,"PcaReport_{dtype}.html"),dtype=dtypes),

        # MultiQC
        join(workpath,"Reports/multiqc_report.html"),

   
        
   rule rawfastqc:
      input: 
        expand("{name}.R1.fastq.gz", name=samples), 
      output: 
        join(workpath,"rawQC")
      priority: 2
      params: 
        rname='pl:rawfastqc',
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        fastqcver=config['bin'][pfamily]['tool_versions']['FASTQCVER'],
      threads: 32
      shell: """
mkdir -p {output};
module load {params.fastqcver};
fastqc {input} -t {threads} -o {output}; 
        """

   rule trim_se:
      input: 
        infq=join(workpath,"{name}.R1."+config['project']['filetype']),
      output: 
        outfq=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
      params: 
        rname='pl:trim_se',
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        cutadaptver=config['bin'][pfamily]['tool_versions']['CUTADAPTVER'],
        fastawithadaptersetd=config['bin'][pfamily]['tool_parameters']['FASTAWITHADAPTERSETD'],
        leadingquality=config['bin'][pfamily]['tool_parameters']['LEADINGQUALITY'],
        trailingquality=config['bin'][pfamily]['tool_parameters']['TRAILINGQUALITY'],
        minlen=config['bin'][pfamily]['tool_parameters']['MINLEN'],
      threads:32
      shell: """
module load {params.cutadaptver};
cutadapt --nextseq-trim=2 --trim-n -n 5 -O 5 -q {params.leadingquality},{params.trailingquality} -m {params.minlen} -b file:{params.fastawithadaptersetd} -j {threads} -o {output.outfq} {input.infq}
"""

   rule fastqc:
      input: 
        expand(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"), name=samples)
      output: 
        join(workpath,"QC")
      priority: 2
      params: 
        rname='pl:fastqc',
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        fastqcver=config['bin'][pfamily]['tool_versions']['FASTQCVER'],
      threads: 32
      shell: """
mkdir -p {output};
module load {params.fastqcver};
fastqc {input} -t {threads} -o {output};
module load python/3.5;
python Scripts/get_read_length.py {output} > {output}/readlength.txt  2> {output}/readlength.err
"""

   rule fastq_screen:
      input: 
        file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
      output: 
        out1=join(workpath,"FQscreen","{name}.R1.trim_screen.txt"),
        out2=join(workpath,"FQscreen","{name}.R1.trim_screen.png"),
        out3=join(workpath,"FQscreen2","{name}.R1.trim_screen.txt"),
        out4=join(workpath,"FQscreen2","{name}.R1.trim_screen.png")
      params: 
        rname='pl:fqscreen',
        batch='--cpus-per-task=24 --mem=64g --time=10:00:00',
        fastq_screen=config['bin'][pfamily]['tool_versions']['FASTQ_SCREEN'],
        outdir = join(workpath,"FQscreen"),
        outdir2 = join(workpath,"FQscreen2"),
        fastq_screen_config=config['bin'][pfamily]['tool_parameters']['FASTQ_SCREEN_CONFIG'],
        fastq_screen_config2=config['bin'][pfamily]['tool_parameters']['FASTQ_SCREEN_CONFIG2'],
        perlver=config['bin'][pfamily]['tool_versions']['PERLVER'],
        bowtie2ver=config['bin'][pfamily]['tool_versions']['BOWTIE2VER'],
      threads: 24
      shell: """
module load {params.bowtie2ver};
module load {params.perlver};
{params.fastq_screen} --conf {params.fastq_screen_config} --outdir {params.outdir} --threads {threads} --subset 1000000 --aligner bowtie2 --force {input.file1}
{params.fastq_screen} --conf {params.fastq_screen_config2} --outdir {params.outdir2} --threads {threads} --subset 1000000 --aligner bowtie2 --force {input.file1}
"""

   rule kraken_se:
      input: 
        fq=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
      output: 
        krakentaxa = join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.taxa.txt"),
        kronahtml = join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.krona.html"),
      params: 
        rname='pl:kraken',
        prefix = "{name}",
        outdir=join(workpath,kraken_dir),
        bacdb=config['bin'][pfamily]['tool_parameters']['KRAKENBACDB'],
        krakenver=config['bin'][pfamily]['tool_versions']['KRAKENVER'],
        kronatoolsver=config['bin'][pfamily]['tool_versions']['KRONATOOLSVER'],
      threads: 24
      shell: """
module load {params.krakenver};
module load {params.kronatoolsver};
cd /lscratch/$SLURM_JOBID;
cp -rv {params.bacdb} /lscratch/$SLURM_JOBID/;
kraken --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` --fastq-input --gzip-compressed --threads {threads} --output /lscratch/$SLURM_JOBID/{params.prefix}.krakenout --preload {input.fq}
kraken-translate --mpa-format --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` /lscratch/$SLURM_JOBID/{params.prefix}.krakenout |cut -f2|sort|uniq -c|sort -k1,1nr > /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa
cut -f2,3 {params.prefix}.krakenout | ktImportTaxonomy - -o {params.prefix}.kronahtml
mv /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa {output.krakentaxa}
mv /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml {output.kronahtml}
"""

   rule star1p:
      input: 
        file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        qcdir=join(workpath,"QC")
        #,file3="fastqc_status_checked.txt"
      output: 
        out1=join(workpath,star_dir,"{name}.SJ.out.tab"),
        out3=temp(join(workpath,star_dir,"{name}.Aligned.out.bam")),
        #out2= "{name}.SJ.out.tab.Pass1.sjdb"
      params: 
        rname='pl:star1p',
        prefix="{name}",
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        starver=config['bin'][pfamily]['tool_versions']['STARVER'],
        stardir=config['references'][pfamily]['STARDIR'],
        filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],
        samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],
        filtertype=config['bin'][pfamily]['FILTERTYPE'],
        filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],
        alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],
        alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],
        filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],
        filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],
        alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],
        alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],
        alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],
        adapter1=config['bin'][pfamily]['ADAPTER1'],
        adapter2=config['bin'][pfamily]['ADAPTER2'],
      threads: 32
      run:
        import glob,json
        rl=int(open(join(input.qcdir,"readlength.txt")).readlines()[0].strip())-1
        dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.stardir+'*/',recursive=False))))
        bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
        cmd="cd {star_dir}; module load "+params.starver+"; STAR --genomeDir "+params.stardir+str(bestdbrl)+" --outFilterIntronMotifs "+params.filterintronmotifs+" --outSAMstrandField "+params.samstrandfield+"  --outFilterType "+params.filtertype+" --outFilterMultimapNmax "+str(params.filtermultimapnmax)+" --alignSJoverhangMin "+str(params.alignsjoverhangmin)+" --alignSJDBoverhangMin "+str(params.alignsjdboverhangmin)+"  --outFilterMismatchNmax "+str(params.filtermismatchnmax)+" --outFilterMismatchNoverLmax "+str(params.filtermismatchnoverlmax)+"  --alignIntronMin "+str(params.alignintronmin)+" --alignIntronMax "+str(params.alignintronmax)+" --alignMatesGapMax "+str(params.alignmatesgapmax)+" --clip3pAdapterSeq "+params.adapter1+" "+params.adapter2+" --readFilesIn "+input.file1+" --readFilesCommand zcat --runThreadN "+str(threads)+" --outFileNamePrefix "+params.prefix+". --outSAMtype BAM Unsorted"
        shell(cmd)
        A=open(join(workpath,"run.json"),'r')
        a=eval(A.read())
        A.close()
        config=dict(a.items())
        config['project']['SJDBOVERHANG']=str(bestdbrl)
        config['project']["STARDIR"]= config['references'][pfamily]['STARDIR']+str(bestdbrl)
        config['project']['READLENGTH']=str(rl+1)
        with open(join(workpath,'run.json'),'w') as F:
          json.dump(config, F, sort_keys = True, indent = 4,ensure_ascii=False)
        F.close()

   rule sjdb:
      input: 
        files=expand(join(workpath,star_dir,"{name}.SJ.out.tab"), name=samples),
      output: 
        out1=join(workpath,star_dir,"uniq.filtered.SJ.out.tab"),
      params: 
        rname='pl:sjdb'
      shell: """
cat {input.files} |sort|uniq|awk -F \"\\t\" '{{if ($5>0 && $6==1) {{print}}}}'|cut -f1-4|sort|uniq|grep \"^chr\"|grep -v \"^chrM\" > {output.out1}
        """

   rule star2p:
      input: 
        file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
        tab=join(workpath,star_dir,"uniq.filtered.SJ.out.tab"),
        qcdir=join(workpath,"QC"),
        #dir="STARINDEX"
      output: 
        out1=temp(join(workpath,star_dir,"{name}.p2.Aligned.sortedByCoord.out.bam")),
        out2=join(workpath,star_dir,"{name}.p2.ReadsPerGene.out.tab"),
        out3=join(workpath,bams_dir,"{name}.p2.Aligned.toTranscriptome.out.bam"),
        out4=join(workpath,star_dir,"{name}.p2.SJ.out.tab"),
        out5=join(workpath,log_dir,"{name}.p2.Log.final.out"),
      params: 
        rname='pl:star2p',
        prefix="{name}.p2",
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        starver=config['bin'][pfamily]['tool_versions']['STARVER'],
        filterintronmotifs=config['bin'][pfamily]['FILTERINTRONMOTIFS'],
        samstrandfield=config['bin'][pfamily]['SAMSTRANDFIELD'],
        filtertype=config['bin'][pfamily]['FILTERTYPE'],
        filtermultimapnmax=config['bin'][pfamily]['FILTERMULTIMAPNMAX'],
        alignsjoverhangmin=config['bin'][pfamily]['ALIGNSJOVERHANGMIN'],
        alignsjdboverhangmin=config['bin'][pfamily]['ALIGNSJDBOVERHANGMIN'],
        filtermismatchnmax=config['bin'][pfamily]['FILTERMISMATCHNMAX'],
        filtermismatchnoverlmax=config['bin'][pfamily]['FILTERMISMATCHNOVERLMAX'],
        alignintronmin=config['bin'][pfamily]['ALIGNINTRONMIN'],
        alignintronmax=config['bin'][pfamily]['ALIGNINTRONMAX'],
        alignmatesgapmax=config['bin'][pfamily]['ALIGNMATESGAPMAX'],
        adapter1=config['bin'][pfamily]['ADAPTER1'],
        adapter2=config['bin'][pfamily]['ADAPTER2'],
        outsamunmapped=config['bin'][pfamily]['OUTSAMUNMAPPED'],
        wigtype=config['bin'][pfamily]['WIGTYPE'],
        wigstrand=config['bin'][pfamily]['WIGSTRAND'], 
        gtffile=config['references'][pfamily]['GTFFILE'], 
        nbjuncs=config['bin'][pfamily]['NBJUNCS'],
        stardir=config['references'][pfamily]['STARDIR']
      threads:32
      run:
        import glob
        rl=int(open(join(input.qcdir,"readlength.txt")).readlines()[0].strip())-1
        dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.stardir+'*/',recursive=False))))
        bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
        cmd="cd {star_dir}; module load "+params.starver+"; STAR --genomeDir "+params.stardir+str(bestdbrl)+" --outFilterIntronMotifs "+params.filterintronmotifs+" --outSAMstrandField "+params.samstrandfield+"  --outFilterType "+params.filtertype+" --outFilterMultimapNmax "+str(params.filtermultimapnmax)+" --alignSJoverhangMin "+str(params.alignsjoverhangmin)+" --alignSJDBoverhangMin "+str(params.alignsjdboverhangmin)+"  --outFilterMismatchNmax "+str(params.filtermismatchnmax)+" --outFilterMismatchNoverLmax "+str(params.filtermismatchnoverlmax)+"  --alignIntronMin "+str(params.alignintronmin)+" --alignIntronMax "+str(params.alignintronmax)+" --alignMatesGapMax "+str(params.alignmatesgapmax)+"  --clip3pAdapterSeq "+params.adapter1+" "+params.adapter2+" --readFilesIn "+input.file1+" --readFilesCommand zcat --runThreadN "+str(threads)+" --outFileNamePrefix "+params.prefix+". --outSAMunmapped "+params.outsamunmapped+" --outWigType "+params.wigtype+" --outWigStrand "+params.wigstrand+" --sjdbFileChrStartEnd "+str(input.tab)+" --sjdbGTFfile "+params.gtffile+" --limitSjdbInsertNsj "+str(params.nbjuncs)+" --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate"
        shell(cmd)
        cmd="sleep 120;cd {workpath};mv {workpath}/{star_dir}/{params.prefix}.Aligned.toTranscriptome.out.bam {workpath}/{bams_dir}; mv {workpath}/{star_dir}/{params.prefix}.Log.final.out {workpath}/{log_dir}"
        shell(cmd)
        

rule samplecondition:
   input: 
    files=expand(join(workpath,star_dir,"{name}.star.count.txt"), name=samples)
   output: 
    out1=join(workpath,star_dir,"sampletable.txt")
   params: 
    rname='pl:samplecondition',
    batch='--mem=4g --time=10:00:00', 
    pathprefix=join(workpath,star_dir),
    groups=config['project']['groups']['rgroups'],
    labels=config['project']['groups']['rlabels'],
    allsamples=config['project']['groups']['rsamps'],
    gtffile=config['references'][pfamily]['GTFFILE']
   run:
        with open(output.out1, "w") as out:
            out.write("sampleName\tfileName\tcondition\tlabel\n")
            i=0
            for f in input.files:
                out.write("{}\t".format(params.allsamples[i]))
                out.write("{}/{}.star.count.txt\t".format(params.pathprefix, params.allsamples[i]))
                out.write("{}\t".format(params.groups[i]))
                out.write("{}\n".format(params.labels[i])) 
                i=i+1
            out.close()

rule get_strandness:
  input: 
    groupsfile=join(workpath,"groups.tab"),
    files=expand(join(workpath,log_dir,"{name}.RnaSeqMetrics.txt"),name=samples),
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
    check_readaccess(input.groupsfile)
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


rule picard:
   input: 
    file1=join(workpath,star_dir,"{name}.p2.Aligned.sortedByCoord.out.bam"),
   output: 
    outstar2=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
    outstar2b=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bai"),
    outstar3=join(workpath,log_dir,"{name}.star.duplic")
   params: 
    rname='pl:picard',
    sampleName="{name}",
    batch='--mem=24g --time=10:00:00 --gres=lscratch:800',
    picardver=config['bin'][pfamily]['tool_versions']['PICARDVER'],
   shell: """
module load {params.picardver};
java -Xmx110g  -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups I={input.file1} O=/lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.bam TMP_DIR=/lscratch/$SLURM_JOBID RGID=id RGLB=library RGPL=illumina RGPU=machine RGSM=sample; 
java -Xmx110g -jar $PICARDJARPATH/picard.jar MarkDuplicates I=/lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.bam O=/lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.dmark.bam TMP_DIR=/lscratch/$SLURM_JOBID CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE={output.outstar3};
mv /lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.dmark.bam {output.outstar2};
mv /lscratch/$SLURM_JOBID/{params.sampleName}.star_rg_added.sorted.dmark.bai {output.outstar2b};
sed -i 's/MarkDuplicates/picard.sam.MarkDuplicates/g' {output.outstar3};
"""

rule preseq:
	params:
		rname = "pl:preseq",
		preseqver=config['bin'][pfamily]['tool_versions']['PRESEQVER'],
	input:
		bam = join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
	output:
		ccurve = join(workpath,preseq_dir,"{name}.ccurve"),
	shell:"""
module load {params.preseqver};
preseq c_curve -B -o {output.ccurve} {input.bam}            
            """


rule stats:
   input: 
    file1=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
   output: 
    outstar1=join(workpath,log_dir,"{name}.RnaSeqMetrics.txt"),
    outstar2=join(workpath,log_dir,"{name}.flagstat.concord.txt"),
   params: 
    rname='pl:stats',
    batch='--mem=24g --time=10:00:00 --gres=lscratch:800',
    picardver=config['bin'][pfamily]['tool_versions']['PICARDVER'],
    samtoolsver=config['bin'][pfamily]['tool_versions']['SAMTOOLSVER'],
    refflat=config['references'][pfamily]['REFFLAT'],
    rrnalist=config['references'][pfamily]['RRNALIST'],
    picardstrand=config['bin'][pfamily]['PICARDSTRAND']
   shell: """
module load R/3.5;
module load {params.picardver};
java -Xmx110g -jar $PICARDJARPATH/picard.jar CollectRnaSeqMetrics REF_FLAT={params.refflat} I={input.file1} O={output.outstar1} RIBOSOMAL_INTERVALS={params.rrnalist}  STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND TMP_DIR=/lscratch/$SLURM_JOBID  VALIDATION_STRINGENCY=SILENT;
sed -i 's/CollectRnaSeqMetrics/picard.analysis.CollectRnaSeqMetrics/g' {output.outstar1}
module load {params.samtoolsver};
samtools flagstat {input.file1} > {output.outstar2};
module load python/3.5;
python Scripts/bam_count_concord_stats.py {input.file1} >> {output.outstar2}
"""


rule prernaseqc:
   input: 
    expand(join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"), name=samples)
   output:
    out1=join(workpath,bams_dir,"files_to_rnaseqc.txt")
   priority: 2
   params: 
    rname='pl:prernaseqc',batch='--mem=4g --time=04:00:00'
   run:
        with open(output.out1, "w") as out:
            out.write("Sample ID\tBam file\tNotes\n")
            for f in input:
                out.write("%s\t"  % f)
                out.write("%s\t"  % f)
                out.write("%s\n"  % f)
            out.close()

rule rnaseqc:
   input:
    join(workpath,bams_dir,"files_to_rnaseqc.txt")
   output:
    join(workpath,"STAR_QC")
   priority: 2
   params: 
    rname='pl:rnaseqc',
    batch='--mem=24g --time=48:00:00',
    bwaver=config['bin'][pfamily]['tool_versions']['BWAVER'],
    rrnalist=config['references'][pfamily]['RRNALIST'],
    rnaseqcver=config['bin'][pfamily]['RNASEQCJAR'],
    gtffile=config['references'][pfamily]['GTFFILE'],
    genomefile=config['references'][pfamily]['GENOMEFILE']
   shell: """
module load {params.bwaver}
var="{params.rrnalist}"
if [  $var == "-" ]; then
      java -Xmx48g -jar {params.rnaseqcver} -n 1000 -s {input} -t {params.gtffile} -r {params.genomefile}  -o {output}
else
      java -Xmx48g -jar {params.rnaseqcver} -n 1000 -s {input} -t {params.gtffile} -r {params.genomefile} -rRNA {params.rrnalist}  -o {output}
fi
"""

rule rnaseq_multiqc:
   input: 
    expand(join(workpath,rseqc_dir,"{name}.Rdist.info"),name=samples),
    expand(join(workpath,"FQscreen","{name}.R1.trim_screen.png"),name=samples),
    expand(join(workpath,log_dir,"{name}.flagstat.concord.txt"),name=samples),
    expand(join(workpath,log_dir,"{name}.RnaSeqMetrics.txt"),name=samples),
    expand(join(workpath,log_dir,"{name}.star.duplic"),name=samples),
    expand(join(workpath,preseq_dir,"{name}.ccurve"),name=samples),
   output: 
    join(workpath,"Reports","multiqc_report.html")
   params: 
    rname="pl:multiqc",
    logsdir=join(workpath,log_dir),
    outdir=join(workpath,"Reports"),
    multiqcver=config['bin'][pfamily]['tool_versions']['MULTIQCVER'],
    qcconfig=config['bin'][pfamily]['CONFMULTIQC']
   threads: 1
   shell: """
module load {params.multiqcver}
cd {params.outdir}
multiqc -f -c {params.qcconfig} --interactive -e cutadapt -d ../
cd {workpath}/slurmfiles
multiqc -f --interactive .
    """

if pe=="yes":

   rule rsem:
      input: 
        file1=join(workpath,bams_dir,"{name}.p2.Aligned.toTranscriptome.out.bam"),
        file2=join(workpath,rseqc_dir,"{name}.strand.info")
      output: 
        out1=join(workpath,degall_dir,"{name}.RSEM.genes.results"),
        out2=join(workpath,degall_dir,"{name}.RSEM.isoforms.results"),
      params:
        rname='pl:rsem',
        prefix="{name}.RSEM",
        outdir=join(workpath,degall_dir),
        batch='--cpus-per-task=16 --mem=32g --time=24:00:00',
        rsemref=config['references'][pfamily]['RSEMREF'],
        rsemver=config['bin'][pfamily]['tool_versions']['RSEMVER'],
        pythonver=config['bin'][pfamily]['tool_versions']['PYTHONVER'],
        annotate=config['references'][pfamily]['ANNOTATE'],
        pythonscript=join(workpath,"Scripts","merge_rsem_results.py"),
      threads: 16
      shell: """
if [ ! -d {params.outdir} ]; then mkdir {params.outdir}; fi
cd {params.outdir}
module load {params.rsemver}
fp=`tail -n1 {input.file2} |awk '{{if($NF > 0.75) print "0.0"; else if ($NF<0.25) print "1.0"; else print "0.5";}}'`
echo $fp
rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  --bam --paired-end -p {threads}  {input.file1} {params.rsemref} {params.prefix} --time --temporary-folder /lscratch/$SLURM_JOBID --keep-intermediate-files --forward-prob=$fp --estimate-rspd
"""


if se=="yes":

   rule rsem:
      input:
        file1=join(workpath,bams_dir,"{name}.p2.Aligned.toTranscriptome.out.bam"),
        file2=join(workpath,rseqc_dir,"{name}.strand.info")
      output:
        out1=join(workpath,degall_dir,"{name}.RSEM.genes.results"),
        out2=join(workpath,degall_dir,"{name}.RSEM.isoforms.results"),
      params: 
        rname='pl:rsem',
        prefix="{name}.RSEM",
        outdir=join(workpath,degall_dir),
        batch='--cpus-per-task=16 --mem=32g --time=24:00:00',
        rsemref=config['references'][pfamily]['RSEMREF'],
        rsemver=config['bin'][pfamily]['tool_versions']['RSEMVER'],
        pythonver=config['bin'][pfamily]['tool_versions']['PYTHONVER'],
        annotate=config['references'][pfamily]['ANNOTATE'],
        pythonscript=join(workpath,"Scripts","merge_rsem_results.py"),
      threads: 16
      shell: """
if [ ! -d {params.outdir} ]; then mkdir {params.outdir}; fi
cd {params.outdir}
module load {params.rsemver}
fp=`tail -n1 {input.file2} |awk '{{if($NF > 0.75) print "0.0"; else if ($NF<0.25) print "1.0"; else print "0.5";}}'`
echo $fp
rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  --bam -p {threads}  {input.file1} {params.rsemref} {params.prefix} --time --temporary-folder /lscratch/$SLURM_JOBID --keep-intermediate-files --forward-prob=$fp --estimate-rspd
"""

rule bam2bw:
	input:
		bam=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
		strandinfo=join(workpath,rseqc_dir,"{name}.strand.info")
	output:
		fbw=join(workpath,bams_dir,"{name}.fwd.bw"),
		rbw=join(workpath,bams_dir,"{name}.rev.bw")
	params:
		rname='pl:bam2bw',
		prefix="{name}",
		deeptoolsver=config['bin'][pfamily]['tool_versions']['DEEPTOOLSVER']
	threads: 32
	shell:"""
module load {params.deeptoolsver};

bam={input.bam}
# Forward strand
bamCoverage -b $bam -o {params.prefix}.fwd.bw --filterRNAstrand forward --binSize 20 --smoothLength 40 -p 32

# Reverse strand
bamCoverage -b $bam -o {params.prefix}.rev.bw --filterRNAstrand reverse --binSize 20 --smoothLength 40 -p 32

# reverse files if method is not dUTP/NSR/NNSR ... ie, R1 in the direction of RNA strand.
fp=`awk '{{if($NF > 0.75) print "0.0"; else if ($NF<0.25) print "D1.0"; else print "0.5";}}' {input.strandinfo}`
if [ $fp -lt 0.25 ];then
mv {params.prefix}.fwd.bw {params.prefix}.fwd.bw.tmp
mv {params.prefix}.rev.bw {params.prefix}.fwd.bw
mv {params.prefix}.fwd.bw.tmp {params.prefix}.rev.bw
fi
"""


rule rsem_merge:
   input:
    files=expand(join(workpath,degall_dir,"{name}.RSEM.genes.results"), name=samples),
    files2=expand(join(workpath,degall_dir,"{name}.RSEM.isoforms.results"), name=samples),
   output: 
    join(workpath,degall_dir,"RSEM.genes.FPKM.all_samples.txt"),
    join(workpath,degall_dir,"RSEM.isoforms.FPKM.all_samples.txt"),
   params: 
    rname='pl:rsem_merge',
    pythonver=config['bin'][pfamily]['tool_versions']['PYTHONVER'],
    annotate=config['references'][pfamily]['ANNOTATE'],
    pythonscript=join(workpath,"Scripts","merge_rsem_results.py"),
   shell: """
module load {params.pythonver}
python {params.pythonscript} {params.annotate} {degall_dir} {degall_dir}
"""


rule rsemcounts:
   input:
    files=expand(join(workpath,degall_dir,"{name}.RSEM.genes.results"), name=samples),
    sampletable=join(workpath,star_dir,"sampletable.txt")
   output: 
    join(workpath,degall_dir,"RawCountFile_RSEM_genes_filtered.txt"),
   params: 
    rname='pl:rsemcounts',
    batch='--mem=8g --time=10:00:00',
    outdir=join(workpath,degall_dir),
    annotate=config['references'][pfamily]['ANNOTATE'],
    rver=config['bin'][pfamily]['tool_versions']['RVER'],
    rscript=join(workpath,"Scripts","rsemcounts.R")
   shell: """
cd {params.outdir}
module load {params.rver}
Rscript {params.rscript} '{params.outdir}' '{input.files}' '{params.annotate}' '{input.sampletable}' 
"""

rule qualicounts:
   input:
    countsmatrix=join(workpath,degall_dir,"RawCountFile_RSEM_genes_filtered.txt"),
    groupsfile=join(workpath,"groups.tab"),
   output:
    outcounts=join(workpath,"QualiMap","RawCountFile_RSEM_genes_filtered_qualimap.txt"),
    globalreport=join(workpath,"QualiMap","GlobalReport.html"),
   params:
    rname='pl:qualicounts',
    sampletable=join(workpath,"QualiMap","qualimap_sample_table.txt"),
    outdir=join(workpath,"QualiMap"),
    info=config['references'][pfamily]['QUALIMAP_INFO'],
   shell: """
module load qualimap/2.2.1
# Remove gene symbols from count matrix
sed 's/|[a-zA-Z0-9]\+//g' {input.countsmatrix} | tail -n +2 > {output.outcounts}
sed '/^$/d' {input.groupsfile} | awk -v OFS='\\t' '{{print $1, $2,"{output.outcounts}", NR+1}}' > {params.sampletable}
qualimap counts -d {params.sampletable} -i {params.info} -outdir {params.outdir}
        """

rule rseqc:
   input: 
    file1=join(workpath,bams_dir,"{name}.star_rg_added.sorted.dmark.bam"),
   output: 
    out1=join(workpath,rseqc_dir,"{name}.strand.info"),
    out4=join(workpath,rseqc_dir,"{name}.Rdist.info")
   params: 
    bedref=config['references'][pfamily]['BEDREF'],
    rseqcver=config['bin'][pfamily]['tool_versions']['RSEQCVER'],
    rname="pl:rseqc"
   shell: """
module load {params.rseqcver}
cd {rseqc_dir}
infer_experiment.py -r {params.bedref} -i {input.file1} -s 1000000 > {output.out1}
read_distribution.py -i {input.file1} -r {params.bedref} > {output.out4}
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
featureCounts -T {threads} -s `cat {input.file2}` -p -t exon -g gene_id -a {params.gtffile} --tmpDir /lscratch/$SLURM_JOBID  -o {output.out}  {input.file1}
sed '1d' {output.out} | cut -f1,7 > {output.res}
"""

rule genecounts: 
   input:
    file1=expand(join(workpath,star_dir,"{name}.star.count.txt"), name=samples),
    file2=join(workpath,star_dir,"sampletable.txt")
   output:
    join(workpath,degall_dir,"RawCountFile_Subread_genes_filtered.txt")
   params: 
    rname='pl:genecounts',
    batch='--mem=8g --time=10:00:00',
    outdir=join(workpath,degall_dir),
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
    join(workpath,degall_dir,"RawCountFile_Subread_junctions_filtered.txt")
   params:
    rname='pl:junctioncounts',
    batch='--mem=8g --time=10:00:00',
    outdir=join(workpath,degall_dir),
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
    join(workpath,degall_dir,"RawCountFile_Subread_genejunctions_filtered.txt")
   params:
    rname='pl:genejunctions',
    batch='--mem=8g --time=10:00:00',
    outdir=join(workpath,degall_dir),
    geneinfo=config['references'][pfamily]['GENEINFO'],
    mincount=config['project']['MINCOUNTGENEJUNCTIONS'],
    minsamples=config['project']['MINSAMPLES'],
    rver=config['bin'][pfamily]['tool_versions']['RVER'],
    rscript=join(workpath,"Scripts","genejunctioncounts.R")
   shell: """
module load {params.rver}
module load bedtools/2.27.1
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
featureCounts -T {threads} -s `cat {input.file2}` -p -t exon -O -g gene_id -a {params.gtffile} --tmpDir /lscratch/$SLURM_JOBID  -o {output.out}  {input.file1}
sed '1d' {output.out} | cut -f1,7 > {output.res}
"""

rule pca:
  input: 
    file1=join(workpath,star_dir,"sampletable.txt"),
    file2=join(workpath,degall_dir,"RawCountFile_{dtype}_filtered.txt"),
  output: 
    outhtml=join(workpath,degall_dir,"PcaReport_{dtype}.html")
  params: 
    rname='pl:pca',
    batch='--mem=24g --time=10:00:00',
    outdir=join(workpath,degall_dir),
    dtype="{dtype}",
    projectId=config['project']['id'],
    projDesc=config['project']['description'].rstrip('\n'),
    rver=config['bin'][pfamily]['tool_versions']['RVER'],
    scripts_dir=join(workpath,"Scripts"),
    rscript1="pcacall.R",
    rscript2="PcaReport.Rmd",
  shell: """
cd {params.outdir}
if [ ! -f {params.rscript1} ]; then cp {params.scripts_dir}/{params.rscript1} {params.outdir}/;fi
if [ ! -f {params.rscript2} ]; then cp {params.scripts_dir}/{params.rscript2} {params.outdir}/;fi
module load {params.rver}
Rscript {params.rscript1} '{params.outdir}' '{output.outhtml}' '{input.file1}' '{input.file2}' '{params.projectId}' '{params.projDesc}' 
"""

