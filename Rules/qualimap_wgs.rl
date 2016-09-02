rule qualimap_wgs:
    input:  "{x}.dedup.bam"
    output: "QC/{x}.qualimapReport","QC/{x}.qualimapReport/genome_results.txt"
    threads: 8
    params: qualimap=config['bin'][pfamily]['QUALIMAP'],organism=config['references'][pfamily]['ORGANISM'],rname="pl:wgs_qualimap"
    shell:  "{params.qualimap} bamqc -bam {input} -c gd {params.organism} -outdir {output} -outformat HTML -nw 500 -p NON-STRAND-SPECIFIC"