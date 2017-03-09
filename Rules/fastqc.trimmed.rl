rule fastqc_trimmed:
    input:  "{x}.R1.trimmed.fastq.gz","{x}.R2.trimmed.fastq.gz"
    output: "QC/{x}.R1.trimmed_fastqc.html","QC/{x}.R2.trimmed_fastqc.html"
    params: fastqc=config['bin'][pfamily]['FASTQC'],adapters=config['references'][pfamily]['fastqc.adapters'],rname="pl:fastqc"
    threads: 8
    shell:  "{params.fastqc} -o QC -f fastq --threads {threads} --contaminants {params.adapters} {input}"
