rule all_initialqc:
    input:  expand("QC/{s}.R1_fastqc.html",s=samples),
#            expand("QC/{s}.R1."+config['project']['filetype']+"_filtered",s=samples),
#            expand("{s}.R1.trimmed.fastq.bz2",s=samples),
#            expand("QC/{s}.qualimapReport",s=samples),
            expand("QC/{s}.{r}.trimmed_fastqc.html",s=samples,r=['R1','R1']),
            expand("QC/{s}.{r}_screen.txt",s=samples,r=['R1','R1']),
            expand("QC/{s}.{r}_screen.png",s=samples,r=['R1','R1']),
            expand("{s}.dedup.bam",s=samples),
            expand("{s}.dedup.bam.onTarget.bam_stats",s=samples),
            expand("{s}.dedup.bam.onTarget.bam",s=samples),
            expand("{s}.dedup.bam.bam_stats",s=samples),
            "multiqc_report.html"
    output: 
