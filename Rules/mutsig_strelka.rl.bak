rule mutsig_strelka:
        input: expand("strelka/oncotator/{p}.maf",p=pairs)
        output: genes="strelka/mutsigCV/somatic.sig_genes.txt",
                maf="strelka/oncotator/strelka_variants.maf"
        params: rname="pl:mutsig_strelka"
        shell: "cat strelka/oncotator/*.maf > {output.maf}; module load MutSig; MutSigCV {output.maf} $MUTSIG_REF/exome_full192.coverage.txt $MUTSIG_REF/gene.covariates.txt strelka/mutsigCV/somatic $MUTSIG_REF/mutation_type_dictionary_file.txt $MUTSIG_REF/chr_files_hg19"