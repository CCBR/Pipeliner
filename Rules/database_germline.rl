rule database_germline:
    input: vcf="exome.relaxedFilter.vcf",
           annotation=config['project']['workpath']+"/full_annot.txt.zip"
    output: dbase="variants.database",
            vcf="exome_genotypes.vcf"
    params: regions=config['references'][pfamily]['REFFLAT'],rname="pl:database"
    shell: "unzip -p full_annot.txt.zip > full_annot.txt; gzip exome.relaxedFilter.vcf; module load samtools; bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' exome.relaxedFilter.vcf.gz > exome_genotypes.vcf; gunzip exome.relaxedFilter.vcf.gz; module load perl/5.18.4; perl Scripts/make_database.pl {input.vcf} {output.vcf}; sed '1!b;s/#//g' variants.database > temp.database; module load R/3.5; Rscript Scripts/dedupcol.R; echo -n "create table vcf (\"" > sqlite3.commands; sed -n 1p temp.database | sed 's/#//g' | sed 's/\t/" text, "/g' | sed '${s/$/\" text);/}' >> sqlite3.commands; echo ".separator \"\\t\"" >> sqlite3.commands; echo ".import temp.database vcf" >> sqlite3.commands; tail -n +2 temp.database > temp.temp.database; mv temp.temp.database temp.database; sqlite3 vcf.db < sqlite3.commands"