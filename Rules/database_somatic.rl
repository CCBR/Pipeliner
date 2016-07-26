rule database_somatic:
    input: vcf="mutect_out/merged_somatic.vcf",
           annotation="full_annot.txt.zip"
    output: dbase="variants.database",
            vcf="exome_genotypes.vcf"
    params: regions=config['references'][pfamily]['REFFLAT'],rname="pl:database"
    shell: "unzip -p full_annot.txt.zip > full_annot.txt; gzip {input.vcf}; module load samtools; bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {input.vcf}.gz > exome_genotypes.vcf; perl Scripts/make_database.pl {input.vcf} {output.vcf}"