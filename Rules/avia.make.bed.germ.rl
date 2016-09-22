rule avia_make_bed_germ:
    input: "combined.vcf"
    output: bed="variants.bed",
            vcf="exome.recode.vcf"
    params: regions="exome_targets.bed",batch ="-l nodes=1:gpfs -q ccr",rname="make_bed"
    shell: """
         module load vcftools; vcftools --vcf combined.vcf --bed {params.regions} --recode --recode-INFO-all --out exome; perl Scripts/avia_make_bed.pl exome.recode.vcf; mkdir sample_vcfs

           """

