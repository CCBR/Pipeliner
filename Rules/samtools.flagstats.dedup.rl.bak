rule flagstats_dedup:
    input:  lambda wildcards: config['project']['units'][wildcards.x]+".dedup.bam"
    output: "QC/{x}.dedup.flags"
    params: sam=config['bin']['SAMTOOLS'],rname="pl:flagstat"
    shell:  "{params.sam} flagstat {input} > {output}"

