rule bwa_mem:
    input:  "{x}.R1.trimmed.fastq.gz","{x}.R2.trimmed.fastq.gz"
    output: "{x}.bam"
    params: bwa=config['bin']['BWA'],ref=config['references']['BWAGENOME'],,rname="pl:bwape",partition="norm", mem="16g", time="24:00:00",sam=config['bin']['SAMTOOLS'],rname="pl:bwamem"
    threads: 4
    shell:  "{params.bwa} mem -t {threads} {params.ref} {input}  > {output}| {params.sam} view -bS - > {output.bam};"

