rule qualimap:
    input:  bam="{x}.dedup.bam",
            bed=config['project']['workpath']+"/exome_targets.bed"
    output: dir="QC/{x}.qualimapReport",txt="QC/{x}.qualimapReport/genome_results.txt"
    threads: 8
    params: qualimap=config['bin'][pfamily]['QUALIMAP'],organism=config['references'][pfamily]['ORGANISM'],regions="exome_targets.bed",rname="pl:qualimap"
    shell:  "module load qualimap/2.2.1; unset DISPLAY; qualimap bamqc -bam {input.bam} --java-mem-size=48G -c gd {params.organism} -ip -outdir {output.dir} -gff {params.regions} -outformat HTML -nt {threads} --skip-duplicated -nw 500 -p NON-STRAND-SPECIFIC"