rule fastqc_fastq:
    input:  "{x}.R1."+config['project']['filetype'],"{x}.R2."+config['project']['filetype']
    output: "QC/{x}.R1_fastqc.html","QC/{x}.R2_fastqc.html"
    params: fastqc=config['bin']['FASTQC'],adapters=config['references']['fastqc.adapters'],rname="pl:fastqc"
    threads: 8
    shell:  "{params.fastqc} -o QC -f fastq --threads {threads} --contaminants {params.adapters} {input}"