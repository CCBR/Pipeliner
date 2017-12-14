num = sum(1 for line in open('pairs'))

if num > 10:
  if config['project']['annotation'] == "hg19":
    rule mutsig_strelka:
        input: config['project']['workpath']+"/strelka_out/oncotator_out/strelka_filtered.maf",
        output: genes=config['project']['workpath']+"/strelka_out/mutsigCV_out/somatic.sig_genes.txt",
        params: rname="pl:mutsig_strelka"
        shell: "module load MutSig; MutSigCV {input} $MUTSIG_REF/exome_full192.coverage.txt $MUTSIG_REF/gene.covariates.txt strelka_out/mutsigCV_out/somatic $MUTSIG_REF/mutation_type_dictionary_file.txt $MUTSIG_REF/chr_files_hg19"

  if config['project']['annotation'] == "hg38":
    rule mutsig_strelka:
        input: config['project']['workpath']+"/strelka_out/oncotator_out/strelka_filtered.maf",
        output: genes=config['project']['workpath']+"/strelka_out/mutsigCV_out/somatic.sig_genes.txt",
        params: rname="pl:mutsig_strelka"
        shell: "echo \"null\" > {output.genes}"

  elif config['project']['annotation'] == "mm10":
    rule mutsig_strelka:
        input: config['project']['workpath']+"/strelka_out/oncotator_out/strelka_filtered.maf",
        output: genes=config['project']['workpath']+"/strelka_out/mutsigCV_out/somatic.sig_genes.txt",
        params: rname="pl:mutsig_strelka"
        shell: "echo \"null\" > {output.genes}"

else:
    rule mutsig_strelka:
        input: config['project']['workpath']+"/strelka_out/oncotator_out/strelka_filtered.maf",
        output: genes=config['project']['workpath']+"/strelka_out/mutsigCV_out/somatic.sig_genes.txt",
        params: rname="pl:mutsig_strelka"
        shell: "echo \"null\" > {output.genes}"