rule avia_make_bed_genomeseq:
    input: "combined.vcf"
    output: bed=config['project']['workpath']+"/variants.bed"
    params: regions=config['references'][pfamily]['REFFLAT'],batch ="-l nodes=1:gpfs -q ccr",rname="make_bed_genome"
    shell: """
         perl Scripts/avia_make_bed.pl {input}; mkdir -p sample_vcfs

           """

