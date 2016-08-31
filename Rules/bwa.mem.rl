rule bwa_mem:
    input:  "{x}.R1.trimmed.fastq.gz","{x}.R2.trimmed.fastq.gz"
    output: bam="{x}.bam"
    params: bwa=config['bin'][pfamily]['BWA'],ref=config['references'][pfamily]['BWAGENOME'],sam=config['bin'][pfamily]['SAMTOOLS'],rname="pl:bwamem"
    threads: 32
    shell:  """
             module load samtools; {params.bwa} mem -t {threads} {params.ref} {input} | samtools view -bS - > {output}

            """
