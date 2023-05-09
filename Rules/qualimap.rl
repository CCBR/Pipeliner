rule qualimap:
    input:  bam="{x}.recal.bam",
            bed=config['project']['workpath']+"/exome_targets.bed"
    output: txt="QC/{x}.qualimapReport/genome_results.txt"
    threads: 8
    params: qualimap=config['bin'][pfamily]['QUALIMAP'],organism=config['references'][pfamily]['ORGANISM'],regions="exome_targets.bed",rname="pl:qualimap", out_dir="QC/{x}.qualimapReport"
    shell:  "module load qualimap/2.2.1; unset DISPLAY; qualimap bamqc -bam {input.bam} --java-mem-size=48G -c gd {params.organism} -ip -outdir {params.out_dir} -gff {params.regions} -outformat HTML -nt {threads} --skip-duplicated -nw 500 -p NON-STRAND-SPECIFIC"
