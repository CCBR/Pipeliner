rule fastqc_rnaseqvar:
    input:  "{x}.R1."+config['project']['filetype'],"{x}.R2."+config['project']['filetype']
    output: "QC/{x}.R1_fastqc.html","QC/{x}.R2_fastqc.html",length="QC/{x}_readlength.txt"
    params: fastqc=config['bin'][pfamily]['FASTQC'],name="{x}",adapters=config['references'][pfamily]['fastqc.adapters'],rname="pl:fastqc"
    threads: 8
    shell:  "{params.fastqc} -o QC -f fastq --threads {threads} --contaminants {params.adapters} {input}; cd QC; perl ../Scripts/get_rnaseqvar_read_length.pl {params.name}"