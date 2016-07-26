rule samtools_flagstats:
    input:  lambda wildcards: config['project']['units'][wildcards.x]+".bam"
    output: "QC/{x}.flags"
    params: sam=config['bin']['SAMTOOLS'],rname="pl:flagstat"
    shell:  "{params.sam} flagstat {input} > {output}"

