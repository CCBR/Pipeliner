from snakemake.utils import R
from os.path import join
configfile: "run.json"
# configfile: "run1.json.1"
from os import listdir
import os

#various subroutines
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

samples=config['project']['groups']['rsamps']

workpath=config['project']['workpath']

# contrastsList = createConstrasts(config['project']['contrasts']['rcontrasts'])
# cpm_cutoff_list = config['project']['contrasts']['rcontrasts_cpm_cutoff']
# mincounts = config['project']['contrasts']['rcontrasts_min_counts']

# contrasts_w_cpm_cutoff_list = []
# for i in zip(contrastsList,cpm_cutoff_list,mincounts):
#   contrasts_w_cpm_cutoff_list.append(str(i[0]) + "_" + str(i[1]) + "_" + str(i[2]))
# print(contrasts_w_cpm_cutoff_list, contrastsList)


#trim #cutadapt
se = ""
pe = ""
workpath = config['project']['workpath']
# 
# if config['project']['nends'] == 2 :
# 	se="no"
# elif config['project']['nends'] == 1 :
# 	se="yes"

se = "yes"
trim_dir='trim'
kraken_dir='kraken'
log_dir="logfiles"
# bams_dir="bams"
# degall_dir="DEG_ALL"
pfamily = config['project']['pfamily']

for d in [trim_dir,kraken_dir,log_dir]:#trim_dir,kraken_dir,bams_dir,star_dir,log_dir,rseqc_dir,preseq_dir,degall_dir
  if not os.path.exists(join(workpath,d)):
    os.mkdir(join(workpath,d))

if not os.path.exists("slurmfiles"):
  os.mkdir ("slurmfiles")
  
if se=="yes":
  if (config['project']['novelMirs'] == "YES"):
    rule all:
      params:
        batch = '--time=168:00:00'
      input:
        join(workpath,'rawQC'),
        join(workpath,'QC'),
        expand(join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.taxa.txt"),name=samples),
        expand(join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.krona.html"),name=samples),
        expand(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),name=samples),
        join(workpath,'Reports/multiqc_report.html'),
        expand(join(workpath,"FQscreen","{name}.R1.trim_screen.txt"),name=samples),
        expand(join(workpath,"FQscreen2","{name}.R1.trim_screen.txt"),name=samples),
        expand(join(workpath,"FQscreen","{name}.R1.trim_screen.png"),name=samples),
        expand(join(workpath,"FQscreen2","{name}.R1.trim_screen.png"),name=samples),
        # # # # expand(join(workpath,"fasta/trim/{name}.trim.fasta.gz"),name=samples),
        # # # # expand(join(workpath,"fasta/trim/{name}.trim.fasta"),name=samples),
        expand(join(workpath,"fasta","{name}.trimfasta.done"),name=samples),
        join(workpath,"fasta","collapsed_trim.fa"),
        join(workpath,"mirdeep2_1p/miRNAs_expressed_all_samples_1p.txt"),
        join(workpath,"fasta/novel_mature.fa"),
        join(workpath,"fasta/novel_hairpin.fa"),
        join(workpath,"mirdeep2_2p/novel_miRNAs_expressed_all_samples_2p.txt"),
        join(workpath,"mirdeep2_results/annotatedCounts.txt"),
        join(workpath,"mirdeep2_results/novelCounts.txt")
  elif(config['project']['novelMirs'] == "NO"):
    rule all:
      params:
        batch = '--time=72:00:00'
      input:
        join(workpath,'rawQC'),
        join(workpath,'QC'),
        expand(join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.taxa.txt"),name=samples),
        expand(join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.krona.html"),name=samples),
        expand(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),name=samples),
        join(workpath,'Reports/multiqc_report.html'),
        expand(join(workpath,"FQscreen","{name}.R1.trim_screen.txt"),name=samples),
        expand(join(workpath,"FQscreen2","{name}.R1.trim_screen.txt"),name=samples),
        expand(join(workpath,"FQscreen","{name}.R1.trim_screen.png"),name=samples),
        expand(join(workpath,"FQscreen2","{name}.R1.trim_screen.png"),name=samples),
        # # # # expand(join(workpath,"fasta/trim/{name}.trim.fasta.gz"),name=samples),
        # # # # expand(join(workpath,"fasta/trim/{name}.trim.fasta"),name=samples),
        expand(join(workpath,"fasta","{name}.trimfasta.done"),name=samples),
        join(workpath,"fasta","collapsed_trim.fa"),
        join(workpath,"mirdeep2_out","miRNAs_expressed_all_samples.txt"),
        join(workpath,"mirdeep2_results/annotated_counts.txt")


  rule raw_fastqc:
    input:
      expand("{name}.R1.fastq.gz", name = samples)
    output:
      join(workpath,"rawQC")
    priority: 2
    params:
      rname = 'pl:rawfastqc',
      batch = '--cpus-per-task=32 --mem=110g --time=48:00:00',
      fastqcver = config['bin'][pfamily]['tool_versions']['FASTQCVER'],
    threads: 32
    shell: """
      mkdir -p {output};
      module load {params.fastqcver};
      fastqc {input} -t {threads} -o {output};
      """
  rule trim_se:
    input:
      infq=join(workpath,"{name}.R1.fastq.gz") #+config['project']['filetype'] used if needed for file flexibility
    output:
      out11=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz")
      #err=join(workpath,"QC","{name}_run_cutadapt.err")
    params:
      rname='pl:trim_se',
      batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
      cutadaptver=config['bin'][pfamily]['tool_versions']['CUTADAPTVER'],
      fastawithadaptersetd=config['bin'][pfamily]['tool_parameters']['FASTAWITHADAPTERSETD'],
      leadingquality=config['bin'][pfamily]['tool_parameters']['LEADINGQUALITY'],
      trailingquality=config['bin'][pfamily]['tool_parameters']['TRAILINGQUALITY'],
      minlen=config['bin'][pfamily]['tool_parameters']['MINLEN']
    threads: 32
    shell:"""
      module load {params.cutadaptver};
      cutadapt --nextseq-trim=2 --trim-n -n 5 -O 5 -q {params.leadingquality},{params.trailingquality} -m {params.minlen} -b file:{params.fastawithadaptersetd} -j {threads} -o {output.out11} {input.infq}
    """
    
  rule kraken:
    input:
      fq=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz")
    output:
      krakentaxa=join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.taxa.txt"),
      kronahtml=join(workpath,kraken_dir,"{name}.trim.fastq.kraken_bacteria.krona.html")
    params:
      rname='pl:kraken',
      prefix="{name}",
      outdir=join(workpath,kraken_dir),
      bacdb=config['bin'][pfamily]['tool_parameters']['KRAKENBACDB'],
      krakenver=config['bin'][pfamily]['tool_versions']['KRAKENVER'],
      kronatoolsver=config['bin'][pfamily]['tool_versions']['KRONATOOLSVER'],
    threads: 24
    shell: """
      module load {params.krakenver};
      module load {params.kronatoolsver};
      if [ ! -d {params.outdir} ]; then mkdir {params.outdir}; fi
      if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID ; fi  
      cd /lscratch/$SLURM_JOBID;
      cp -rv {params.bacdb} /lscratch/$SLURM_JOBID/;
      dbName=`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'`;
      kraken --db /lscratch/$SLURM_JOBID/$dbName --fastq-input --gzip-compressed --threads {threads} --output /lscratch/$SLURM_JOBID/{params.prefix}.krakenout --preload {input.fq}
      kraken-translate --mpa-format --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` /lscratch/$SLURM_JOBID/{params.prefix}.krakenout |cut -f2|sort|uniq -c|sort -k1,1nr > /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa
      cut -f2,3 /lscratch/$SLURM_JOBID/{params.prefix}.krakenout | ktImportTaxonomy - -o /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml
      mv /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa {output.krakentaxa}
      mv /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml {output.kronahtml}
    """
      
  rule fastqc:
    input:
      expand(join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),name=samples),
    output:
      join(workpath,"QC")
    priority: 2
    params:
      rname='pl:fastqc',
      batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
      fastqcver=config['bin'][pfamily]['tool_versions']['FASTQCVER'],
    threads: 32
    shell: """
      if [ ! -d {output} ]; then mkdir -p {output}; fi
      module load {params.fastqcver};
      fastqc {input} -t {threads} -o {output};
      module load python/3.5;
      python Scripts/get_read_length.py {output} > {output}/readlength.txt 2> {output}/readlength.err
      """
    
  rule fastq_screen:
    input:
      file1=join(workpath,trim_dir,"{name}.R1.trim.fastq.gz")
    output:
      out1=join(workpath,"FQscreen","{name}.R1.trim_screen.txt"),
      out2=join(workpath,"FQscreen","{name}.R1.trim_screen.png"),
      out3=join(workpath,"FQscreen2","{name}.R1.trim_screen.txt"),
      out4=join(workpath,"FQscreen2","{name}.R1.trim_screen.png"),
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
  rule multiqc:
    input:
      expand(join(workpath,"FQscreen","{name}.R1.trim_screen.png"),name=samples),
    output:
      join(workpath,"Reports","multiqc_report.html")
    params:
      rname="pl:multiqc",
      logsdir=join(workpath,log_dir),
      outdir=join(workpath,"Reports"),
      multiqcver=config['bin'][pfamily]['tool_versions']['MULTIQCVER'],
      qcconfig=config['bin'][pfamily]['tool_parameters']['CONFMULTIQC']
    threads: 1
    shell: """
      module load {params.multiqcver}
      cd {params.outdir}
      multiqc -f -c {params.qcconfig} --interactive -e cutadapt -d ../
      cd {workpath}/slurmfiles
      multiqc -f --interactive .
    """

  rule fastaConversion:
    input:
      join(workpath,trim_dir,"{name}.R1.trim.fastq.gz"),
    output:
      join(workpath,"fasta","{name}.trimfasta.done")
    params:
      rname='pl:fastq2fasta',
      batch='--cpus-per-task=32 --mem=64g --time=24:00:00',
      fastxtoolkitver = config['bin'][pfamily]['tool_versions']['FASTXTOOLKITVER'],
      outfile = join(workpath,"fasta/trim/{name}.trim.fasta"),
      outdir = join(workpath,"fasta","trim")

    threads: 32
    shell: """
      module load {params.fastxtoolkitver};
      mkdir -p {params.outdir}
      zcat {input} | fastq_to_fasta -n -v -o {params.outfile}
      sed -i -e 's/\s/|/g' {params.outfile}
      touch {output}
    """

  rule mirdeep2_mapper:
    input:
      configFile=join(workpath,"mirdeep_config.txt"),
      # file1=expand(join(workpath,"fasta/trim/{name}.trim.fasta.gz"),name=samples)
      fileCatch=expand(join(workpath,"fasta","{name}.trimfasta.done"),name=samples)
      # fileCatch=join(workpath,"fasta/trim/{name}.trimfasta.done")
    output:
      #join(workpath,"fasta/trim/{name}.trim.fasta.gz"),
      collapsedFile=join(workpath,"fasta","collapsed_trim.fa"),
      arfFile=join(workpath,"trimmed.arf")      
    params:
      rname='pl:mirdeep2_mapper',
      #outdir = join(workpath,"fasta/collapsed"),
      mirdeepver=config['bin'][pfamily]['tool_versions']['MIRDEEP2VER'],
      bowtieindex = config['references'][pfamily]['BOWTIE_REF'],
      trimmedFolder = join(workpath,"fasta","trim"),
      out_tar = join(workpath,"fasta","trim.tar.gz"),
      # file1=expand(join(workpath,"fasta/trim/{name}.trim.fasta"),name=samples)   
    threads: 8
    shell: """
      module load {params.mirdeepver}
      mapper.pl {input.configFile} -d -c -j -l 18 -m -q -p {params.bowtieindex} -s {output.collapsedFile} -t {output.arfFile} -v -o 8 2>> logfiles/mapper.log
      tar -czvf {params.out_tar} {params.trimmedFolder} --remove-files
    """
    
  rule mirdeep2_1p:
    input:
      collapsedFile = join(workpath,"fasta","collapsed_trim.fa"),
      arfFile = join(workpath,"trimmed.arf")
    output:
      annotAlign=join(workpath,"mirdeep2_1p","miRNAs_expressed_all_samples_1p.txt"),
      novelmaturefasta = join(workpath,"fasta/novel_mature.fa"),
      novelhairpinfasta = join(workpath,"fasta/novel_hairpin.fa")
    params:
      rname="pl:mirdeep2_1p",
      mirdeepver=config['bin'][pfamily]['tool_versions']['MIRDEEP2VER'],
      genomefasta = config['references'][pfamily]['FASTA_REF'],
      maturefasta = config['references'][pfamily]['MIRBASE_MATURE'],
      hairpinfasta = config['references'][pfamily]['MIRBASE_HAIRPIN'],
      organism=config['references'][pfamily]['ORGANISM'],
      resultdir=join(workpath,"mirdeep2_1p"),
        
    threads: 2
    shell: """
      module load {params.mirdeepver}
      mkdir -p /lscratch/$SLURM_JOBID
      mkdir -p {params.resultdir}
      cd {params.resultdir}
      sed -e 's/\s.*$//g' {params.genomefasta} > /lscratch/$SLURM_JOBID/genome.fa #remove whitespace, per miRDeep2
      sed -e 's/Homo sapiens/Homo_sapiens/g' {params.maturefasta}| sed -e 's/\s/|/g' > /lscratch/$SLURM_JOBID/mature.fa
      sed -e 's/Homo sapiens/Homo_sapiens/g' {params.hairpinfasta}| sed -e 's/\s/|/g' > /lscratch/$SLURM_JOBID/hairpin.fa
      
      miRDeep2.pl {input.collapsedFile} /lscratch/$SLURM_JOBID/genome.fa {input.arfFile} /lscratch/$SLURM_JOBID/mature.fa none /lscratch/$SLURM_JOBID/hairpin.fa -t {params.organism} -P
      cat mirna_results*/novel_mature_*.fa <(sed '/^>.*/s/$/_star/' mirdeep2_1p/mirna_results*/novel_star*.fa) > {output.novelmaturefasta}
      cp mirna_results*/novel_pres_*.fa {output.novelhairpinfasta}
      mv miRNAs_expressed_all*.csv {output.annotAlign}
    """
  rule mirdeep2_2p:
    input:
      collapsedFile = join(workpath,"fasta","collapsed_trim.fa"),
      arfFile = join(workpath,"trimmed.arf"),
      novelmaturefasta = join(workpath,"fasta","novel_mature.fa"),
      novelhairpinfasta = join(workpath,"fasta","novel_hairpin.fa")
    output:
      join(workpath,"mirdeep2_2p","novel_miRNAs_expressed_all_samples_2p.txt")
    params:
      rname="pl:mirdeep2_2p",
      mirdeepver=config['bin'][pfamily]['tool_versions']['MIRDEEP2VER'],
      genomefasta = config['references'][pfamily]['FASTA_REF'],
      organism=config['references'][pfamily]['ORGANISM'],
      resultdir=join(workpath,"mirdeep2_2p")   
    threads: 16
    shell: """
      module load {params.mirdeepver}
      mkdir -p /lscratch/$SLURM_JOBID
      mkdir -p {params.resultdir}
      cd {params.resultdir}
      sed -e 's/\s.*$//g' {params.genomefasta} > /lscratch/$SLURM_JOBID/genome.fa #remove whitespace, per miRDeep2
      quantifier.pl -r {input.collapsedFile} -m {input.novelmaturefasta} -p {input.novelhairpinfasta} -T 8 
      mv miRNAs_expressed_all*.csv {output}
    """
  rule capture_counts_2p:
    input:
      mirbaseResults=join(workpath,"mirdeep2_1p/miRNAs_expressed_all_samples_1p.txt"),
      novelResults=join(workpath,"mirdeep2_2p/novel_miRNAs_expressed_all_samples_2p.txt")
    output:
      mirbaseCts=join(workpath,"mirdeep2_results/annotatedCounts.txt"),
      novelCts=join(workpath,"mirdeep2_results/novelCounts.txt")
    params:
      rname="pl:capture_counts_2p",
      configFile=join(workpath,"mirdeep_config.txt"),
      resultdir = join(workpath,"mirdeep2_results"),
      sampleCounts = len(samples)+4
    threads: 4
    shell: """
      mkdir -p {params.resultdir}
      cut -f 1-{params.sampleCounts} {input.mirbaseResults} > {output.mirbaseCts}
      cut -f 1-{params.sampleCounts} {input.novelResults} > {output.novelCts}
    """
  #   
  rule mirdeep2_single:
    input:
      collapsedFile = join(workpath,"fasta/collapsed_trim.fa"),
      configFile=join(workpath,"mirdeep_config.txt"),
    output:
      join(workpath,"mirdeep2_out","miRNAs_expressed_all_samples.txt")
    params:
      rname="pl:mirdeep2_single",
      maturefasta = config['references'][pfamily]['MIRBASE_MATURE'],
      hairpinfasta = config['references'][pfamily]['MIRBASE_HAIRPIN'],
      mirdeepver=config['bin'][pfamily]['tool_versions']['MIRDEEP2VER'],
      genomefasta = config['references'][pfamily]['FASTA_REF'],
      organism=config['references'][pfamily]['ORGANISM'],
      resultdir=join(workpath,"mirdeep2_out")
    threads:16
    shell: """
      module load {params.mirdeepver}
      mkdir -p /lscratch/$SLURM_JOBID
      mkdir -p {params.resultdir}
      cd {params.resultdir}
      #sed -e 's/\s.*$//g' {params.genomefasta} > /lscratch/$SLURM_JOBID/genome.fa #remove whitespace, per miRDeep2
      sed -e 's/Homo sapiens/Homo_sapiens/g' {params.maturefasta}| sed -e 's/\s/|/g' > /lscratch/$SLURM_JOBID/mature.fa
      sed -e 's/Homo sapiens/Homo_sapiens/g' {params.hairpinfasta}| sed -e 's/\s/|/g' > /lscratch/$SLURM_JOBID/hairpin.fa
      quantifier.pl -c {input.configFile} -r {input.collapsedFile} -m /lscratch/$SLURM_JOBID/mature.fa -p /lscratch/$SLURM_JOBID/hairpin.fa -T 8 -t {params.organism} -P
      mv miRNAs_expressed_all*.csv {output}
    """
  
  rule capture_counts_single:
    input:
      mirbaseResults=join(workpath,"mirdeep2_out/miRNAs_expressed_all_samples.txt"),
    output:
      mirbaseCts=join(workpath,"mirdeep2_results/annotated_counts.txt"),
    params:
      rname="pl:capture_counts_single",
      configFile=join(workpath,"mirdeep_config.txt"),
      resultdir = join(workpath,"mirdeep2_results"),
      sampleCounts = len(samples)+4,
      groupsFile = join(workpath,"groups.tab")
    threads: 4
    shell: """
      mkdir -p {params.resultdir}
      cut -f 1-{params.sampleCounts} {input.mirbaseResults} > {output.mirbaseCts}
      while read FILE ABB; do FILE=${{FILE%%.trim*}};FILE=${{FILE##*trim/}}; sed -i "s/$ABB/$FILE/" {output.mirbaseCts}; done < {params.configFile}
      while read FILE GROUP NAME; do sed -i "s/$FILE/$NAME/" {output.mirbaseCts}; done < {params.groupsFile}
      sed -i 's/^#//' {output.mirbaseCts}
    """
    
  # 
    

#needed rules:
#QC Side
#trim: cutadapt - done
#kraken: -done
#
#post-trim fastqc- done
#multiqc - done
#Alignment
#fastq to fasta:fastxtoolkit - done

#collapse (miRDeep2) - done
#align to mature miRs (Bowtie1) -done
#align to precursor miRs (Bowtie1) - done
#align to genome (Bowtie1/miRDeep2 pass1) - done
#align to novel (Bowtie1/miRDeep2 pass2) - done
#expand counts (miRDeep2)
#ShortStack integration
#Differential expression analysis (edgeR)
