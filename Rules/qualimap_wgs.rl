rule qualimap_wgs:
    input:  "{x}.dedup.bam"
    output: "QC/{x}.qualimapReport","QC/{x}.qualimapReport/genome_results.txt"
    threads: 8
    params: qualimap=config['bin'][pfamily]['QUALIMAP'],organism=config['references'][pfamily]['ORGANISM'],rname="pl:wgs_qualimap"
    shell:  "module load qualimap/2.2.1; unset DISPLAY; qualimap bamqc -nt $SLURM_CPUS_PER_TASK --java-mem-size=96G -bam {input} -c gd {params.organism} -outdir {output} -ip --skip-duplicated -nt {threads} -outformat HTML -nw 500 -p NON-STRAND-SPECIFIC"