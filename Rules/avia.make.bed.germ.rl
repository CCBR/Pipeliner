rule avia_make_bed_germ:
    input: "combined.vcf"
    output: bed="variants.bed",
            vcf="exome.recode.vcf"
    params: regions="exome_targets.bed",batch ="-l nodes=1:gpfs -q ccr",rname="make_bed"
    threads: 8
    shell: """
         module load GATK/3.5-0; GATK -m 24G SelectVariants -V combined.vcf -o {output.vcf} -L {params.regions} -R /fdb/GATK_resource_bundle/b37/human_g1k_v37.fasta -nt 8; perl Scripts/avia_make_bed.pl exome.recode.vcf; mkdir -p sample_vcfs

           """

