rule fastq_screen:
    input:  "{x}.R1."+config['project']['filetype'],"{x}.R2."+config['project']['filetype']
    output: "QC/{x}.R1_screen.txt",
            "QC/{x}.R1_screen.png",
            "QC/{x}.R2_screen.txt",
            "QC/{x}.R2_screen.png"
    params: fastq_screen=config['bin'][pfamily]['FASTQ_SCREEN'],
            outdir = "QC",
            config=config['references'][pfamily]['FASTQ_SCREEN_CONFIG'],rname="pl:fastq_screen"
    threads: 8
    shell:  "module load bowtie; {params.fastq_screen} --conf {params.config} --outdir {params.outdir} --subset 0 --aligner bowtie2 --force {input}"