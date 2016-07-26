rule bwa_mem:
    input:  "{x}.R1.trimmed.fastq.gz","{x}.R2.trimmed.fastq.gz"
    output: bam="{x}.sam"
    params: bwa=config['bin'][pfamily]['BWA'],ref=config['references'][pfamily]['BWAGENOME'],partition="ccr", mem="128g", time="48:00:00",sam=config['bin'][pfamily]['SAMTOOLS'],rname="pl:bwamem"
    threads: 8
    shell:  """
             {params.bwa} mem -t {threads} {params.ref} {input}  > {output}

            """
