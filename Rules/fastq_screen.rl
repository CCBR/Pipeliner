rule fastq_screen:
    input:  "{x}.R1.trimmed.fastq.gz","{x}.R2.trimmed.fastq.gz"
    output: "QC/{x}.R1.trimmed_screen.txt",
            "QC/{x}.R1.trimmed_screen.png",
            "QC/{x}.R2.trimmed_screen.txt",
            "QC/{x}.R2.trimmed_screen.png"
    params: fastq_screen=config['bin'][pfamily]['FASTQ_SCREEN'],
            outdir = "QC",
            config=config['references'][pfamily]['FASTQ_SCREEN_CONFIG'],rname="pl:fastq_screen"
    threads: 8
    shell:  "module load bowtie; {params.fastq_screen} --conf {params.config} --outdir {params.outdir} --subset 1000000 --aligner bowtie2 --force {input}"