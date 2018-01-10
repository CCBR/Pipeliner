rule multiqc:
    input: expand("QC/{s}.qualimapReport",s=samples), expand("{s}.dedup.bam.onTarget.bam_stats",s=samples), expand("{s}.dedup.bam.bam_stats",s=samples), expand("QC/{s}_run_trimmomatic.err",s=samples), expand("QC/{s}.qualimapReport/genome_results.txt",s=samples), expand("{s}.sorted.txt",s=samples), expand("QC/{s}.R1.trimmed_screen.png", s=samples), expand("QC/{s}.R2.trimmed_screen.png", s=samples)
    output: "multiqc_report.html"
    params: fastqc=config['bin'][pfamily]['FASTQC'],adapters=config['references'][pfamily]['fastqc.adapters'],rname="pl:multiqc"
    threads: 1
    shell:  "module load multiqc/1.3; multiqc -f ."


