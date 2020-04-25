rule avia_make_bed_germ:
    input: "exome.relaxedFilter.vcf"
    output: bed=config['project']['workpath']+"/variants.bed",
    params: regions="exome_targets.bed",genome=config['references'][pfamily]['GENOME'],batch ="-l nodes=1:gpfs -q ccr",rname="make_bed"
    shell: """
         module load perl/5.18.4; perl Scripts/avia_make_bed.pl {input}; mkdir -p sample_vcfs

           """

