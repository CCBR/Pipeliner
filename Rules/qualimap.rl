rule qualimap:
    input:  bam="{x}.dedup.bam",
            bed=config['project']['workpath']+"/exome_targets.bed"
    output: "QC/{x}.qualimapReport","QC/{x}.qualimapReport/genome_results.txt"
    threads: 8
    params: qualimap=config['bin'][pfamily]['QUALIMAP'],organism=config['references'][pfamily]['ORGANISM'],regions="exome_targets.bed",rname="pl:qualimap"
#    shell:  "{params.qualimap} bamqc -bam {input} -c gd {params.organism} -outfile {output} -gff {params.regions} -outformat HTML -nw 500 -p NON-STRAND-SPECIFIC -nt {threads}"
    shell:  "module load qualimap/2.2.1; qualimap bamqc --java-mem-size=96G -bam {input.bam} -c gd {params.organism} -ip -outdir {output} -gff {params.regions} -outformat HTML -nt {threads} --skip-duplicated -nw 500 -p NON-STRAND-SPECIFIC"