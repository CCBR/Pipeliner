rule all-rnaseqvar-qc:
    input:  expand("QC/{s}.R1_fastqc.html",s=samples),
            expand("QC/{s}.{r}.trimmed_fastqc.html",s=samples,r=['R1','R1']),
            expand("QC/{s}.{r}.trimmed_screen.txt",s=samples,r=['R1','R1']),
            expand("QC/{s}.{r}.trimmed_screen.png",s=samples,r=['R1','R1']),
            expand("{s}.dedup.bam",s=samples),
#            expand("{s}.dedup.bam.onTarget.bam_stats",s=samples),
#            expand("{s}.dedup.bam.onTarget.bam",s=samples),
#            expand("{s}.dedup.bam.bam_stats",s=samples),
            "multiqc_report.html",
            expand("QC/{s}_run_trimmomatic.err",s=samples),
            expand("QC/{s}.qualimapReport/genome_results.txt",s=samples),
#            config['project']['id']+"_"+config['project']['flowcellid']+".xlsx"
    output:
    params: rname="final"
    shell:  """
             mv *.out slurmfiles/; module load perl/5.18.4; perl Scripts/summarize_usage.pl; rm *sorted.bam.bai; mv *sorted.txt logfiles/

            """