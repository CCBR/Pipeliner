rule avia_make_bed_germ:
    input: "combined.vcf"
    output: bed="variants.bed",
            vcf="exome.recode.vcf"
    params: regions="exome_targets.bed",genome=config['references'][pfamily]['GENOME'],batch ="-l nodes=1:gpfs -q ccr",rname="make_bed"
    shell: """
         module load GATK/3.5-0; GATK -m 24G SelectVariants -V combined.vcf -o {output.vcf} -L {params.regions} -R {params.genome}; perl Scripts/avia_make_bed.pl exome.recode.vcf; mkdir -p sample_vcfs

           """

