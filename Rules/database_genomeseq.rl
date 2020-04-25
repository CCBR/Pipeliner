rule database_genomeseq:
    input: vcf="combined.relaxedFilter.vcf",
           annotation=config['project']['workpath']+"/full_annot.txt.zip"
    output: dbase="variants.database",
            vcf="combined_genotypes.vcf"
    params: regions=config['references'][pfamily]['REFFLAT'],rname="pl:database"
    shell: "unzip -p {input.annotation} > full_annot.txt; gzip {input.vcf}; module load samtools; bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' combined.vcf.gz > {output.vcf}; gunzip combined.vcf.gz; module load perl/5.18.4; perl Scripts/make_database.pl {input.vcf} {output.vcf}"