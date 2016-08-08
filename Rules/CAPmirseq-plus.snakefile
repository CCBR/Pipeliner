import os

configfile: "run.json"

#l=list(config['sample_info']['units'].values())
#l=list(config['project']['samples'].values())
#samples=list(config['units'].values())
#samples=[i[0].split('.fastq')[0] for i in l]

SAMPLES=":".join(config['project']['contrasts']['rsamps'])
GROUPS=":".join(config['project']['contrasts']['rgroups'])
PAIRS=":".join(config['project']['contrasts']['rcontrasts'])
samples=config['project']['contrasts']['rsamps']

rule mirseq_final:
    input: expand("{out}/init.done",out=config['project']['workpath']),
           expand("{out}/fastqs/{x}.cutadapt.fastq",x=samples,out=config['project']['workpath']),
           expand("{out}/qc/fastqc_pretrim/{x}/{x}_fastqc.zip",x=samples,out=config['project']['workpath']),
           expand("{out}/qc/fastqc_posttrim/{x}/{x}.cutadapt_fastqc.zip",x=samples,out=config['project']['workpath']),
            expand("{out}/bams-bwa/{x}.bam",x=samples,out=config['project']['workpath']),
            expand("{out}/bams/{x}.bam",x=samples,out=config['project']['workpath']),
            expand("{out}/mirdeep2/{x}/{x}.reads.fa",x=samples,out=config['project']['workpath']),
            expand("{out}/mirdeep2/{x}/mirdeep2.log",x=samples,out=config['project']['workpath']),
            config['project']['workpath']+"/variants/mirna_variants.vcf",
            expand("{p}/qc/other_rna/{x}.gencode.genecount.txt",x=samples,p=config['project']['workpath']),
            config['project']['workpath']+"/expression/mature_miRNA_expression.xls",
            config['project']['workpath']+"/differential_expression/expression_boxplots.pdf",
            config['project']['workpath']+"/SampleSummary.xls",
            config['project']['workpath']+"/MainDocument.html",
            expand("{out}/mirspring/{x}.html",x=samples,out=config['project']['workpath']),
            expand("{out}/mirspring/{x}_bwa.html",x=samples,out=config['project']['workpath']),            

rule mirseq_init:
    input:
    output: expand("{p}/init.done",p=config['project']['workpath'])
    params: out=config['project']['workpath'],mem="16G",time="4:00:00",partition="ccr",rname="mir:init"
    priority: 50
    shell: """
	mkdir -p {params.out}
	mkdir -p {params.out}/logs
	mkdir -p {params.out}/bams
	mkdir -p {params.out}/igv
	mkdir -p {params.out}/mirdeep2
	mkdir -p {params.out}/expression
	mkdir -p {params.out}/variants
	mkdir -p {params.out}/qc
	mkdir -p {params.out}/qc/fastqc_pretrim
	mkdir -p {params.out}/qc/fastqc_posttrim
	mkdir -p {params.out}/qc/other_rna
	mkdir -p {params.out}/differential_expression
	mkdir -p {params.out}/fastqs
	mkdir -p {params.out}/config
        touch {params.out}/init.done
           """

## rule trimmomatic:
##     input:  config['project']['workpath']+"/{x}.fastq.gz"
##     output: config['references'][pfamily]['OUTPUT_DIR']+"/fastqs/{x}.cutadapt.fastq"
##     params: trimmomatic=config['references'][pfamily]['tool_info']['tool_parameters']['TRIMMOMATIC'],
##             adapterfile=config['references'][pfamily]['tool_info']['tool_parameters']['TRIMMOMATIC.ADAPTERS'],mem="16G",time="4:00:00",partition="ccr",name="mir:trimmomatic"
##     threads: 4
##     shell:  """
##             {params.trimmomatic} SE -threads {threads} -phred33 {input} {output} ILLUMINACLIP:{params.adapterfile}:3:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:17
##            """
## 
rule mirseq_cutadapt:
       input: config['project']['workpath']+"/{x}.fastq.gz"
       output: config['project']['workpath']+"/fastqs/{x}.cutadapt.fastq",config['project']['workpath']+"/fastqs/{x}.tooshort.fastq",config['project']['workpath']+"/fastqs/{x}.cutadapt.log",
       params: cpath=config['bin'][pfamily]['tool_paths']['CUTADAPT_PATH'],cparams=config['bin'][pfamily]['tool_parameters']['CUTADAPT_PARAMS'],qtrim='20',indir=config['project']['workpath'],out=config['project']['workpath']+"/fastqs",mem="32G",time="4:00:00",partition="ccr",rname="mir:cutadapt"
       threads: 1
       shell: """
   	{params.cpath}/cutadapt {params.cparams} -q {params.qtrim} {input} -o {output[0]} --too-short-output={output[1]} > {output[2]}
             """ 

rule mirseq_fastqc_pretrim:
    input: config['project']['workpath']+"/{x}.fastq.gz"
    output: config['project']['workpath']+"/qc/fastqc_pretrim/{x}/{x}_fastqc.zip"
    params: fastqc=config['bin'][pfamily]['tool_paths']['FASTQC_PATH'],indir=config['project']['workpath'],out=config['project']['workpath']+"/qc/fastqc_pretrim",mem="16G",time="4:00:00",partition="ccr",rname="mir:pretrim"
    threads: 1    
    shell: """
	mkdir -p {params.out}/{wildcards.x}
	{params.fastqc}/fastqc -o {params.out}/{wildcards.x} {params.indir}/{wildcards.x}.fastq.gz
           """

rule mirseq_fastqc_posttrim:
    input: config['project']['workpath']+"/fastqs/{x}.cutadapt.fastq"
    output: config['project']['workpath']+"/qc/fastqc_posttrim/{x}/{x}.cutadapt_fastqc.zip"
    params: fastqc=config['bin'][pfamily]['tool_paths']['FASTQC_PATH'],indir=config['project']['workpath']+"/fastqs",out=config['project']['workpath']+"/qc/fastqc_posttrim",mem="16G",time="4:00:00",partition="ccr",rname="mir:posttrim"
    threads: 1    
    shell: """
	mkdir -p {params.out}/{wildcards.x}
	{params.fastqc}/fastqc -o {params.out}/{wildcards.x} {params.indir}/{wildcards.x}.cutadapt.fastq
           """


rule mirseq_mirdeep2_mapper:
    input: config['project']['workpath']+"/fastqs/{x}.cutadapt.fastq"
    output: config['project']['workpath']+"/mirdeep2/{x}/{x}.reads.fa",config['project']['workpath']+"/mirdeep2/{x}/{x}.reads_vs_genome.arf"
    params: mirdeep=config['bin'][pfamily]['tool_paths']['MIRDEEP2_PATH'],indir=config['project']['workpath']+"/fastqs",out=config['project']['workpath']+"/mirdeep2",mapper_params=config['bin'][pfamily]['tool_parameters']['MAPPER_PARAMS'],bowtie_ref=config['references'][pfamily]['reference_files']['BOWTIE_REF'],mem="16G",time="4:00:00",partition="ccr",rname="mir:mapper"
    threads: 1
    shell: """
        export PATH=/data/dwheeler/CAP-miRSEQ/bin:$PATH
	mkdir -p {params.out}/{wildcards.x}
	{params.mirdeep}/mapper.pl {input} {params.mapper_params} -p {params.bowtie_ref} -s {params.out}/{wildcards.x}/{wildcards.x}.reads.fa -t {params.out}/{wildcards.x}/{wildcards.x}.reads_vs_genome.arf -o {threads}
           """

rule mirseq_make_bams:
    input: config['project']['workpath']+"/fastqs/{x}.cutadapt.fastq"
    output: config['project']['workpath']+"/bams/{x}.bam",config['project']['workpath']+"/bams/{x}.bowtie.log"
    params: out=config['project']['workpath']+"/bams",bowtie=config['bin'][pfamily]['tool_paths']['BOWTIE_PATH'],bowtie_ref=config['references'][pfamily]['reference_files']['BOWTIE_REF'],bowtie_params=config['bin'][pfamily]['tool_parameters']['BOWTIE_PARAMS'],quals="--phred33-quals",samtools=config['bin'][pfamily]['tool_paths']['SAMTOOLS_PATH'],addorreplacereadgroups_params=config['bin'][pfamily]['tool_parameters']['ADDORREPLACEREADGROUPS_PARAMS'],java_path=config['bin'][pfamily]['tool_paths']['JAVA_PATH'],picard_path=config['bin'][pfamily]['tool_paths']['PICARD_PATH'],addorreplacereadgroups_jvm_mem=config['bin'][pfamily]['java_parameters']['ADDORREPLACEREADGROUPS_JVM_MEM'],script_path=config['bin'][pfamily]['tool_paths']['SCRIPT_PATH'],mem="32G",time="4:00:00",partition="ccr",rname="mir:bowtie"
    threads: 1    
    shell: """
	{params.bowtie}/bowtie {params.quals} {params.bowtie_params} --sam-RG ID:{wildcards.x} --sam-RG SM:{wildcards.x} {params.bowtie_ref} {input} {params.out}/{wildcards.x}.aligned.sam 2> {params.out}/{wildcards.x}.bowtie.log

	{params.samtools}/samtools view -SH {params.out}/{wildcards.x}.aligned.sam > {params.out}/{wildcards.x}.aligned.mq.sam

	{params.samtools}/samtools view -S {params.out}/{wildcards.x}.aligned.sam | {params.script_path}/dw_setqual.pl >> {params.out}/{wildcards.x}.aligned.mq.sam
	
	rm {params.out}/{wildcards.x}.aligned.sam

	{params.java_path}/java {params.addorreplacereadgroups_jvm_mem} -jar {params.picard_path}/AddOrReplaceReadGroups.jar INPUT={params.out}/{wildcards.x}.aligned.mq.sam OUTPUT={params.out}/{wildcards.x}.bam SORT_ORDER=coordinate TMP_DIR={params.out} RGID={wildcards.x} RGPU={wildcards.x} RGSM={wildcards.x} {params.addorreplacereadgroups_params}

	rm {params.out}/{wildcards.x}.aligned.mq.sam
	{params.samtools}/samtools index {params.out}/{wildcards.x}.bam
        {params.samtools}/samtools flagstat {params.out}/{wildcards.x}.bam >{params.out}/{wildcards.x}.bam.flagstat
           """

rule mirseq_make_bams_bwa:
#    input: config['project']['workpath']+"/{x}.fastq.gz"
    input: config['project']['workpath']+"/fastqs/{x}.cutadapt.fastq"
    output: config['project']['workpath']+"/bams-bwa/{x}.bam"
    params: out=config['project']['workpath']+"/bams-bwa",samtools=config['bin'][pfamily]['tool_paths']['SAMTOOLS_PATH'],genome=config['references'][pfamily]['reference_files']['REF_GENOME'],mem="32G",time="4:00:00",partition="ccr",rname="mir:bwa"
    threads: 4    
    shell: """
    /usr/local/apps/bwa/0.7.10/bwa mem -t {threads} {params.genome} {input} > {params.out}/{wildcards.x}.sam

     {params.samtools}/samtools view -Shu {params.out}/{wildcards.x}.sam > {output}
     {params.samtools}/samtools flagstat {params.out}/{wildcards.x}.bam > {params.out}/{wildcards.x}.bam.flagstat

           """



rule mirseq_mirdeep2:
    input: config['project']['workpath']+"/mirdeep2/{x}/{x}.reads.fa",config['project']['workpath']+"/mirdeep2/{x}/{x}.reads_vs_genome.arf"
    output: config['project']['workpath']+"/mirdeep2/{x}/mirdeep2.log"
    params: script_path=config['bin'][pfamily]['tool_paths']['SCRIPT_PATH'],mature=config['references'][pfamily]['reference_files']['MIRBASE_MATURE'],precursor=config['references'][pfamily]['reference_files']['MIRBASE_HAIRPIN'],mirdeep2_params=config['bin'][pfamily]['tool_parameters']['MIRDEEP2_PARAMS'],mirdeep2_close_species=config['bin'][pfamily]['tool_parameters']['MIRDEEP2_CLOSE_SPECIES'],genome=config['references'][pfamily]['reference_files']['REF_GENOME'],out=config['project']['workpath']+"/mirdeep2",mirdeep2_path=config['bin'][pfamily]['tool_paths']['MIRDEEP2_PATH'],bowtie_path=config['bin'][pfamily]['tool_paths']['BOWTIE_PATH'],squid_path=config['bin'][pfamily]['tool_paths']['SQUID_PATH'],vienna_path=config['bin'][pfamily]['tool_paths']['VIENNA_PATH'],randfold_path=config['bin'][pfamily]['tool_paths']['RANDFOLD_PATH'],mem="16G",time="8:00:00",partition="ccr",rname="mir:mdeep"
    threads: 1    
    shell: """

module load viennarna
module load bowtie/1.1.1
module load randfold
module load mirdeep
export PATH=/data/dwheeler/CAP-miRSEQ/bin:$PATH
cd {params.out}/{wildcards.x} && {params.script_path}/miRDeep2.pl {wildcards.x}.reads.fa {params.genome} {wildcards.x}.reads_vs_genome.arf {params.mature} {params.mirdeep2_close_species} {params.precursor} {params.mirdeep2_params} 2>&1|tee mirdeep2.log

#	rm -Rf ./dir_prepare_signature*
#	rm -Rf ./expression_analyses
#	rm -Rf ./mirdeep_runs
#	rm $collapsed_seqs_fa
#	rm $mapped_reads_arf

           """

rule mirseq_variants:
    input: expand("{p}/bams/{x}.bam",x=samples,p=config['project']['workpath'])
    output: config['project']['workpath']+"/variants/mirna_variants.vcf"
    params: genome=config['references'][pfamily]['reference_files']['REF_GENOME_IUPAC'],out=config['project']['workpath']+"/variants",script_path=config['bin'][pfamily]['tool_paths']['SCRIPT_PATH'],mirgff=config['references'][pfamily]['reference_files']['MIRBASE_GFF'],bedpath=config['bin'][pfamily]['tool_paths']['BEDTOOLS_PATH'],gatkjar=config['bin'][pfamily]['tool_paths']['GATK_JAR'],unifiedgenotyper_params=config['bin'][pfamily]['tool_parameters']['UNIFIEDGENOTYPER_PARAMS'],vcftools_path=config['bin'][pfamily]['tool_paths']['VCFTOOLS_PATH'],units=expand("{s}",s=samples),java_path=config['bin'][pfamily]['tool_paths']['JAVA_PATH'],unifiedgenotyper_jvm_mem=config['bin'][pfamily]['java_parameters']['UNIFIEDGENOTYPER_JVM_MEM'],vcftools_perllib=config['bin'][pfamily]['tool_paths']['VCFTOOLS_PERLLIB'],mem="16G",time="4:00:00",partition="ccr",rname="mir:variants"
    threads: 1    
    run: 
        I=" -I "+" -I ".join(input)
        cmd="""
        {params.script_path}/dw_makebeds.sh {params.mirgff} {params.out} {params.bedpath}
        {params.java_path}/java {params.unifiedgenotyper_jvm_mem} -jar {params.gatkjar} -T UnifiedGenotyper -R {params.genome} {I} -L {params.out}/mirbase_precursor.tmp.bed -o {params.out}/mirna_variants.vcf {params.unifiedgenotyper_params}
        {params.script_path}/dw_annotation_report.sh {params.out} {params.bedpath} {params.vcftools_path} {params.vcftools_perllib}

        {params.script_path}/dw_write_header.sh {params.out} {params.units}
        {params.script_path}/dw_merge_annotations.sh {params.out}

#	rm $output_dir/mirbase_precursor.tmp.bed
#	rm $output_dir/mirbase_mature.tmp.bed
#	rm $output_dir/seed_region.tmp.bed
#	rm $output_dir/mirna_variants.vcf.txt.pre
#	rm $output_dir/mirna_variants.vcf.txt.pre.mat
#	rm $output_dir/mirna_variants.vcf.txt.pre.mat.seed
#	rm $output_dir/mirna_variants.vcf.txt.bed.precursor.bed
#	rm $output_dir/mirna_variants.vcf.txt.bed.mature.bed
#	rm $output_dir/mirna_variants.vcf.txt.bed.seed.bed
#	rm $output_dir/mirna_variants.vcf.txt
#	rm $output_dir/mirna_variants.vcf.txt.bed
#
"""
        shell(cmd)     

rule mirseq_gencode_classification:
    input: config['project']['workpath']+"/bams/{x}.bam"
    output: config['project']['workpath']+"/qc/other_rna/{x}.gencode.genecount.txt"
    params: java_path=config['bin'][pfamily]['tool_paths']['JAVA_PATH'],sortsam_jvm_mem=config['bin'][pfamily]['java_parameters']['SORTSAM_JVM_MEM'],picard_path=config['bin'][pfamily]['tool_paths']['PICARD_PATH'],input_dir=config['project']['workpath']+"/bams",out=config['project']['workpath']+"/qc/other_rna",sortsam_params=config['bin'][pfamily]['tool_parameters']['SORTSAM_PARAMS'],htseq_path=config['bin'][pfamily]['tool_paths']['HTSEQ_PATH'],htseq_params=config['bin'][pfamily]['tool_parameters']['HTSEQ_PARAMS'],gencode_gtf=config['references'][pfamily]['reference_files']['GENCODE_GTF'],script_path=config['bin'][pfamily]['tool_paths']['SCRIPT_PATH'],python_path=config['bin'][pfamily]['tool_paths']['PYTHON_PATH'],htseq_lib_path=config['bin'][pfamily]['tool_paths']['HTSEQ_LIB_PATH'],samtools_path=config['bin'][pfamily]['tool_paths']['SAMTOOLS_PATH'],rscript_path=config['bin'][pfamily]['tool_paths']['RSCRIPT_PATH'],mem="16G",time="4:00:00",partition="ccr",rname="mir:classify"
    threads: 1    
    shell: """

{params.java_path}/java {params.sortsam_jvm_mem} -jar {params.picard_path}/SortSam.jar INPUT={params.input_dir}/{wildcards.x}.bam OUTPUT={params.out}/{wildcards.x}.queryname.bam SORT_ORDER=queryname TMP_DIR={params.out}/ {params.sortsam_params}

{params.samtools_path}/samtools view {params.out}/{wildcards.x}.queryname.bam | {params.htseq_path}/htseq-count {params.htseq_params} - {params.gencode_gtf} > {params.out}/{wildcards.x}.gencode.genecount.txt

{params.script_path}/dw_gencode.sh {params.python_path} {params.htseq_lib_path} {params.gencode_gtf} {params.out} {wildcards.x} {params.script_path} {params.rscript_path}

           """

rule mirseq_expression_reports:
    input: expand("{out}/mirdeep2/{x}/mirdeep2.log",x=samples,out=config['project']['workpath'])
    output: config['project']['workpath']+"/expression/mature_miRNA_expression.xls",config['project']['workpath']+"/expression/miRNA_expression_raw.xls"

    params: script_path=config['bin'][pfamily]['tool_paths']['SCRIPT_PATH'],input_dir=config['project']['workpath']+"/mirdeep2",out=config['project']['workpath']+"/expression",samples=SAMPLES,mem="16G",time="4:00:00",partition="ccr",rname="mir:expression"
    threads: 1    
    run: 
       shell("{params.script_path}/dw_expression_reports.sh {params.input_dir} {params.out} {params.samples} {params.script_path}")

        
rule mirseq_differential_expression:
    input: config['project']['workpath']+"/expression/mature_miRNA_expression.xls"
    output: config['project']['workpath']+"/differential_expression/expression_boxplots.pdf"
    params: script_path=config['bin'][pfamily]['tool_paths']['SCRIPT_PATH'],out=config['project']['workpath']+"/differential_expression",diff_expression=config['references'][pfamily]['run_info']['DIFF_EXPRESSION'],diff_expression_analyses=config['references'][pfamily]['run_info']['DIFF_EXPRESSION_ANALYSES'],samples=SAMPLES,groups=GROUPS,pairs=PAIRS,mem="16G",time="4:00:00",partition="ccr",rname="mir:diffexp"
    threads: 1    
    shell: """
#       units=":".join(samples)
{params.script_path}/dw_differential_expression.sh  {input} {params.out} {params.diff_expression} {params.diff_expression_analyses} {params.script_path} {params.samples} {params.groups} {params.pairs}


           """
rule mirseq_summarize:
    input: expand("{p}/fastqs/{x}.cutadapt.log",x=samples,p=config['project']['workpath']),expand("{p}/bams/{x}.bowtie.log",x=samples,p=config['project']['workpath']),expand("{p}/expression/miRNA_expression_raw.xls",p=config['project']['workpath']),
#expand("{p}/mirdeep2/{x}/expression_*.html",x=samples,p=config['project']['workpath']),expand("{p}/{x}.precursor.reads.txt",x=samples,p=config['project']['workpath'])
    output: config['project']['workpath']+"/SampleSummary.xls"
    params: out=config['project']['workpath'],trim_adapter=config['references'][pfamily]['run_info']['TRIM_ADAPTER'],script_path=config['bin'][pfamily]['tool_paths']['SCRIPT_PATH'],samples=SAMPLES,mem="16G",time="4:00:00",partition="ccr",rname="mir:summarize"
    threads: 1    
    run: 
       units=":".join(samples)
       shell("{params.script_path}/dw_sample_summary.sh {params.out} {params.samples} {params.trim_adapter}")

           
rule mirseq_main_document:
    input: config['project']['workpath']+"/SampleSummary.xls"
    output: config['project']['workpath']+"/MainDocument.html"
    params: out=config['project']['workpath'],trim_adapter=config['references'][pfamily]['run_info']['TRIM_ADAPTER'],flowcell=config['references'][pfamily]['run_info']['FLOWCELL'],tool=config['references'][pfamily]['run_info']['TOOL'],call_snvs=config['references'][pfamily]['run_info']['CALL_SNVS'],diff_expression=config['references'][pfamily]['run_info']['DIFF_EXPRESSION'],diff_expression_analyses=config['references'][pfamily]['run_info']['DIFF_EXPRESSION_ANALYSES'],email=config['references'][pfamily]['run_info']['EMAIL'],script_path=config['bin'][pfamily]['tool_paths']['SCRIPT_PATH'],delivery_folder=config['project']['workpath'],tool_info=config['references'][pfamily]['run_info']['TOOL_INFO'],genome_build=config['references'][pfamily]['run_info']['GENOME_BUILD'],server="",samples=SAMPLES,mem="4G",time="4:00:00",partition="ccr",rname="mir:main"
    threads: 1    
    run: 
       units=":".join(samples)

       if (config['references'][pfamily]['run_info']['CALL_SNVS'] == "YES" ):
           snvs=1
       if (config['references'][pfamily]['run_info']['TRIM_ADAPTER'] == "YES" ):
           trim=1
       if (config['references'][pfamily]['run_info']['DIFF_EXPRESSION'] == "YES" ):
           diff=config['references'][pfamily]['run_info']['DIFF_EXPRESSION_ANALYSES']
       O=open(config['project']['workpath']+"/pfamily.tmp","w")
       I=eval(open("run.json","r").read())
       for k in I['references'][pfamily]['run_info'].keys():
           O.write("{0}={1}\n".format(k,I['references'][pfamily]['run_info'][k]))
       O.close()

       shell("{params.script_path}/dw_main_document.sh {params.out} {params.script_path} {params.flowcell} {params.tool} {params.call_snvs} {params.trim_adapter} {params.diff_expression} {params.diff_expression_analyses} {params.email};perl {params.script_path}/dw_create_igv.pl {params.out}/igv {params.samples} {params.delivery_folder} {params.tool_info} {params.server} {params.genome_build};cp {params.script_path}/IGV_Setup.doc {params.out}/igv;perl {params.script_path}/dw_main_document.pl {params.out}/pfamily.tmp {params.out}/MainDocument.html {params.out}/SampleSummary.xls {snvs} {trim} {diff};cp {params.script_path}/CAP-miRSeq_workflow.png {params.out}")


rule mirseq_miRspring_bowtie:
    input: config['project']['workpath']+"/bams/{x}.bam"
    output: config['project']['workpath']+"/mirspring/{x}.html"
    params: script_path=config['bin'][pfamily]['tool_paths']['SCRIPT_PATH'],output=config['project']['workpath'],mirbase_files=config['references'][pfamily]['reference_files']['MIRBASE_FILES'],mem="16G",time="4:00:00",partition="ccr",rname="mir:mirspring"
    shell: """
           if [ ! -d {params.output}/mirspring ]
           then
           mkdir {params.output}/mirspring
           fi

           cd {params.script_path} && perl {params.script_path}/BAM_to_Intermediate.pl -ml 0 -s hsa -pre {params.mirbase_files}/hsa.35nt.fasta -gff {params.mirbase_files}/hsa_dw.gff3 -mat {params.mirbase_files}/mature.fa -bam {params.output}/bams/{wildcards.x}.bam -ref 0 -out {params.output}/mirspring/{wildcards.x}.tmp.txt -flank 35

           cd {params.script_path} && perl {params.script_path}/Intermediate_to_miRspring.pl -in {params.output}/mirspring/{wildcards.x}.tmp.txt -s hsa -ref 0 -mat {params.mirbase_files}/mature.fa -out {params.output}/mirspring/{wildcards.x}.html -flank 35
           """

rule mirseq_miRspring_bwa:
    input: config['project']['workpath']+"/bams-bwa/{x}.bam"
    output: config['project']['workpath']+"/mirspring/{x}_bwa.html"
    params: script_path=config['bin'][pfamily]['tool_paths']['SCRIPT_PATH'],output=config['project']['workpath'],mirbase_files=config['references'][pfamily]['reference_files']['MIRBASE_FILES'],mem="16G",time="4:00:00",partition="ccr",rname="mir:mirspring"
    shell: """
           if [ ! -d {params.output}/mirspring ]
           then
           mkdir {params.output}/mirspring
           fi

           cd {params.script_path} && perl {params.script_path}/BAM_to_Intermediate.pl -ml 0 -s hsa -pre {params.mirbase_files}/hsa.35nt.fasta -gff {params.mirbase_files}/hsa_dw.gff3 -mat {params.mirbase_files}/mature.fa -bam {params.output}/bams-bwa/{wildcards.x}.bam -ref 0 -out {params.output}/mirspring/{wildcards.x}_bwa.tmp.txt -flank 35

           cd {params.script_path} && perl {params.script_path}/Intermediate_to_miRspring.pl -in {params.output}/mirspring/{wildcards.x}_bwa.tmp.txt -s hsa -ref 0 -mat {params.mirbase_files}/mature.fa -out {params.output}/mirspring/{wildcards.x}_bwa.html -flank 35
           """
rule mirseq_targetscan:
    input:
    output:
    params: script_path=config['bin'][pfamily]['tool_paths']['SCRIPT_PATH'],output=config['project']['workpath'],target_data=config['references'][pfamily]['reference_files']['TARGET_DATA'],mem="4G",time="4:00:00",partition="ccr",rname="mir:tscan"
    shell: """

           {params.script_path}/targetscan_70.pl {params.target_data}/miR_Family_info_sample.txt {params.target_data}/UTR_sequences_all.txt {params.output}/targetscan_70_output.txt


           """
