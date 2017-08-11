rule avia_make_bed_rnaseqvar:
    input: "combined.strictFilter.vcf"
    output: bed=config['project']['workpath']+"/variants.bed"
    params: genome=config['references'][pfamily]['GENOME'],batch ="-l nodes=1:gpfs -q ccr",rname="make_bed"
    shell: """
         perl Scripts/avia_make_bed.pl {input}

           """

