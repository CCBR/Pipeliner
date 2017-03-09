rule avia_make_bed_somatic:
    input: "mutect2/merged_somatic.vcf"
    output:"variants.bed"
    params: batch ="-l nodes=1:gpfs -q ccr",rname="avia.make.bed.somatic"
    shell: """
         perl Scripts/avia_make_bed.pl {input}

           """

