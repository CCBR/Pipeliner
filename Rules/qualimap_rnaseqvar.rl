rule qualimap_rna:
    input:  bam="{x}.dedup.bam"
    output: "QC/{x}.qualimapReport","QC/{x}.qualimapReport/genome_results.txt"
    threads: 8
    params: qualimap=config['bin'][pfamily]['QUALIMAP'],organism=config['references'][pfamily]['ORGANISM'],regions=config['references'][pfamily]['GTFFILE'],rname="pl:qualimap"
    shell:  "{params.qualimap} bamqc -bam {input.bam} -c gd {params.organism} -outdir {output} -gff {params.regions} -outformat HTML -nw 500 -p NON-STRAND-SPECIFIC"