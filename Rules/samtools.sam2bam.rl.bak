rule samtools_sam2bam:
    input:  "{x}.sam",
    output: temp("{x}.bam")
    params: sam=config['bin']['SAMTOOLS'],rname="pl:sam2bam"
    shell:  "{params.sam} view -bS {input} > {output}"

