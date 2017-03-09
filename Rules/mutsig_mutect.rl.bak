rule mutsig_mutect:
        input: expand("mutect/oncotator/{p}.maf",p=pairs)
        output: genes="mutect/mutsigCV/somatic.sig_genes.txt",
                maf="mutect/oncotator/mutect_variants.maf"
        params: rname="pl:mutsig_mutect"
        shell: "cat mutect/oncotator/*.maf > {output.maf}; module load MutSig; MutSigCV {output.maf} $MUTSIG_REF/exome_full192.coverage.txt $MUTSIG_REF/gene.covariates.txt mutect/mutsigCV/somatic $MUTSIG_REF/mutation_type_dictionary_file.txt $MUTSIG_REF/chr_files_hg19"