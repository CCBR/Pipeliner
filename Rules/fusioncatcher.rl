rule starfusion:
     input:  file1="{x}.R1.trimmed.fastq.gz",file2="{x}.R2.trimmed.fastq.gz",file3="QC/{x}.R1.trimmed_fastqc.html",file4="QC/{x}.R2.trimmed_fastqc.html",length="QC/{x}_readlength.txt"
     output: out1=temp("{x}Aligned.out.bam"),out4="{x}SJ.out.tab",out5="{x}Log.progress.out"
     params: rname='starfusion',sample="{x}"
     threads: 32
     shell: module load STAR/2.5.2b; mkdir starfusion/{params.sample}; 