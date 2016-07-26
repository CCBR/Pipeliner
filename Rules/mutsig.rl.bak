rule mutsig:
        input: expand("mutect2/oncotator/{p}.maf",p=pairs)
        output: genes="mutect2/mutsigCV/somatic.sig_genes.txt",
                maf="mutect2/oncotator/mutect2_variants.maf"
        params: rname="pl:mutsig"
        shell: "cat mutect2/oncotator/*.maf > {output.maf}; module load MutSig; MutSigCV {output.maf} $MUTSIG_REF/exome_full192.coverage.txt $MUTSIG_REF/gene.covariates.txt mutect2/mutsigCV/somatic $MUTSIG_REF/mutation_type_dictionary_file.txt $MUTSIG_REF/chr_files_hg19"