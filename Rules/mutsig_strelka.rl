rule mutsig_strelka:
        input: expand("strelka_out/oncotator_out/{p}.maf",p=pairs)
        output: genes="strelka_out/mutsigCV_out/somatic.sig_genes.txt",
                maf="strelka_out/oncotator_out/strelka_variants.maf"
        params: rname="pl:mutsig_strelka"
        shell: "cat strelka_out/oncotator_out/*.maf > {output.maf}; module load MutSig; MutSigCV {output.maf} $MUTSIG_REF/exome_full192.coverage.txt $MUTSIG_REF/gene.covariates.txt strelka_out/mutsigCV_out/somatic $MUTSIG_REF/mutation_type_dictionary_file.txt $MUTSIG_REF/chr_files_hg19"