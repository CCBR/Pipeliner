rule avia_make_bed_somatic:
    input: config['project']['workpath']+"/mutect2_out/merged_somatic.vcf"
    output:config['project']['workpath']+"/variants.bed"
    params: batch ="-l nodes=1:gpfs -q ccr",rname="avia.make.bed.somatic"
    shell: """
         module load perl/5.18.4; perl Scripts/avia_make_bed.pl {input}

           """

