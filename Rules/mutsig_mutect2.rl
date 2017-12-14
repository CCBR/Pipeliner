num = sum(1 for line in open('pairs'))

if num > 10:
  if config['project']['annotation'] == "hg19":
    rule mutsig_mutect2:
        input: config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_filtered.maf",
        output: genes=config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
        params: rname="pl:mutsig_mutect2"
        shell: "module load MutSig; MutSigCV {input} $MUTSIG_REF/exome_full192.coverage.txt $MUTSIG_REF/gene.covariates.txt mutect2_out/mutsigCV_out/somatic $MUTSIG_REF/mutation_type_dictionary_file.txt $MUTSIG_REF/chr_files_hg19"
  
  elif config['project']['annotation'] == "hg38":
    rule mutsig_mutect2:
        input: config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_filtered.maf",
        output: genes=config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
        params: rname="pl:mutsig_mutect2"
        shell: "echo \"null\" > {output.genes}"

  elif config['project']['annotation'] == "mm10":
    rule mutsig_mutect2:
        input: config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_filtered.maf",
        output: genes=config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
        params: rname="pl:mutsig_mutect2"
        shell: "echo \"null\" > {output.genes}"

else:
    rule mutsig_mutect2:
        input: config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_filtered.maf",
        output: genes=config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
        params: rname="pl:mutsig_mutect2"
        shell: "echo \"null\" > {output.genes}"