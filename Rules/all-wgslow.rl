rule all_wgslow:
    input:  expand("QC/{s}.R1_fastqc.html",s=samples),
            expand("QC/ngsqc/{s}.R1."+config['project']['filetype']+"_filtered",s=samples),
            expand("QC/{s}.qualimapReport",s=samples),
            expand("QC/{s}.{r}.trimmed_fastqc.html",s=samples,r=['R1','R1']),
            expand("{s}.dedup.bam",s=samples),
