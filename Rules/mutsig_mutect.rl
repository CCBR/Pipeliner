rule mutsig_mutect:
        input: expand("mutect_out/oncotator_out/{p}.maf",p=pairs)
        output: genes="mutect_out/mutsigCV_out/somatic.sig_genes.txt",
                maf="mutect_out/oncotator_out/mutect_variants.maf"
        params: rname="pl:mutsig_mutect"
        shell: "cat mutect_out/oncotator_out/*.maf > {output.maf}; module load MutSig; MutSigCV {output.maf} $MUTSIG_REF/exome_full192.coverage.txt $MUTSIG_REF/gene.covariates.txt mutect_out/mutsigCV_out/somatic $MUTSIG_REF/mutation_type_dictionary_file.txt $MUTSIG_REF/chr_files_hg19"