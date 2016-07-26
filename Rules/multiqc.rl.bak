rule multiqc:
    input: expand("QC/{s}.qualimapReport",s=samples)
    output: "multiqc_report.html"
    params: fastqc=config['bin']['FASTQC'],adapters=config['references']['fastqc.adapters'],rname="pl:multiqc"
    threads: 1
    shell:  "module load multiqc; multiqc -f -e featureCounts ."


