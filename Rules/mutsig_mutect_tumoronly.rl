import glob
samps = len(glob.glob1(".","*.recal.bam"))

if samps>10:
  if config['project']['annotation'] == "hg19":
    rule mutsig_mutect_tumoronly:
        input: config['project']['workpath']+"/mutect_out/oncotator_out/final_filtered.maf",
        output: genes=config['project']['workpath']+"/mutect_out/mutsigCV_out/somatic.sig_genes.txt",
        params: rname="pl:mutsig"
        shell: "module load MutSig; MutSigCV {input} $MUTSIG_REF/exome_full192.coverage.txt $MUTSIG_REF/gene.covariates.txt mutect_out/mutsigCV_out/somatic $MUTSIG_REF/mutation_type_dictionary_file.txt $MUTSIG_REF/chr_files_hg19"
  
  elif config['project']['annotation'] == "hg38":
    rule mutsig_mutect_tumoronly:
        input: config['project']['workpath']+"/mutect_out/oncotator_out/final_filtered.maf",
        output: genes=config['project']['workpath']+"/mutect_out/mutsigCV_out/somatic.sig_genes.txt",
        params: rname="pl:mutsig"
        shell: "echo \"null\" > {output.genes}"

  elif config['project']['annotation'] == "mm10":
    rule mutsig_mutect_tumoronly:
        input: config['project']['workpath']+"/mutect_out/oncotator_out/final_filtered.maf",
        output: genes=config['project']['workpath']+"/mutect_out/mutsigCV_out/somatic.sig_genes.txt",
        params: rname="pl:mutsig"
        shell: "echo \"null\" > {output.genes}"

else:
    rule mutsig_mutect_tumoronly:
        input: config['project']['workpath']+"/mutect_out/oncotator_out/final_filtered.maf",
        output: genes=config['project']['workpath']+"/mutect_out/mutsigCV_out/somatic.sig_genes.txt",
        params: rname="pl:mutsig"
        shell: "echo \"null\" > {output.genes}"