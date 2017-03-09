rule database_somatic:
    input: vcf=config['project']['workpath']+"/mutect_out/merged_somatic.vcf",
           annotation=config['project']['workpath']+"/full_annot.txt.zip"
    output: dbase=config['project']['workpath']+"/variants.database",
            vcf=config['project']['workpath']+"/exome_genotypes.vcf"
    params: regions=config['references'][pfamily]['REFFLAT'],rname="pl:database"
    shell: "unzip -p full_annot.txt.zip > full_annot.txt; gzip {input.vcf}; module load samtools; bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' {input.vcf}.gz > exome_genotypes.vcf; gunzip {input.vcf}.gz; perl Scripts/make_database.pl {input.vcf} {output.vcf}"