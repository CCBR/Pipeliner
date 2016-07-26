rule avia_make_bed_germ:
    input: "combined.vcf"
    output: bed="variants.bed",
            vcf="exome.recode.vcf"
    params: regions=config['references'][pfamily]['REFFLAT'],batch ="-l nodes=1:gpfs -q ccr",rname="make_bed"
    shell: """
         module load vcftools; vcftools --vcf combined.vcf --bed {params.regions} --recode --recode-INFO-all --out exome; perl Scripts/avia_make_bed.pl exome.recode.vcf

           """

