if config['project']['annotation'] == "hg19":
  rule oncotator_mutect:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           vcf=config['project']['workpath']+"/mutect_out/{x}.FINAL.vcf",
           dir=config['project']['workpath']+"/mutect_out"
    output: maf=config['project']['workpath']+"/mutect_out/oncotator_out/{x}.maf"
    params: rname="pl:oncotator_mutect"
    shell: "module load oncotator; oncotator -v -o TCGAMAF -i VCF -c $ONCOTATOR_DATASOURCE/tx_exact_uniprot_matches.txt -a Tumor_Sample_Barcode:{input.tumor} -a Matched_Norm_Sample_Barcode:{input.normal} --skip-no-alt --db-dir /fdb/oncotator/oncotator_v1_ds_Jan262015 {input.vcf} {output.maf} hg19"

else:
  rule oncotator_mutect:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           vcf=config['project']['workpath']+"/mutect_out/{x}.FINAL.vcf",
           dir=config['project']['workpath']+"/mutect_out"
    output: maf=config['project']['workpath']+"/mutect_out/oncotator_out/{x}.maf"
    params: rname="pl:oncotator_mutect"
    shell: "echo \"null\" > {output.maf}"