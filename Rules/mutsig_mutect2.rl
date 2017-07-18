if config['project']['annotation'] == "hg19":
  rule mutsig_mutect2:
        input: expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{p}.maf",p=pairs)
        output: genes=config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
                maf=config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_variants.maf"
        params: rname="pl:mutsig_mutect2"
        shell: "cat mutect2_out/oncotator_out/*.maf > {output.maf}; module load MutSig; MutSigCV {output.maf} $MUTSIG_REF/exome_full192.coverage.txt $MUTSIG_REF/gene.covariates.txt mutect2_out/mutsigCV_out/somatic $MUTSIG_REF/mutation_type_dictionary_file.txt $MUTSIG_REF/chr_files_hg19"
  
else:

  rule mutsig_mutect2:
        input: expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{p}.maf",p=pairs)
        output: genes=config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
                maf=config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_variants.maf"
        params: rname="pl:mutsig_mutect2"
        shell: "echo \"null\" > {output.maf}; echo \"null\" > {output.genes}"