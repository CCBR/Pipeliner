rule maftools:
  input: mutect2=config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_variants.maf",
         mutect=config['project']['workpath']+"/mutect_out/oncotator_out/mutect_variants.maf",
         strelka=config['project']['workpath']+"/strelka_out/oncotator_out/strelka_variants.maf",
  output: mutect2summary=config['project']['workpath']+"/mutect2_out/mutect2_maf_summary.pdf",
          mutect2oncoprint=config['project']['workpath']+"/mutect2_out/mutect2_oncoplot.pdf",
          mutectsummary=config['project']['workpath']+"/mutect_out/mutect_maf_summary.pdf",
          mutectoncoprint=config['project']['workpath']+"/mutect_out/mutect_oncoplot.pdf",
          strelkasummary=config['project']['workpath']+"/strelka_out/strelka_maf_summary.pdf",
          strelkaoncoprint=config['project']['workpath']+"/strelka_out/strelka_oncoplot.pdf",
  params: dir=config['project']['workpath'],rname="pl:maftools"
  shell: "perl Scripts/prep_mafs.pl {input.mutect2} mutect2; perl Scripts/prep_mafs.pl {input.mutect} mutect; perl Scripts/prep_mafs.pl {input.strelka} strelka; module load R/3.4.0_gcc-6.2.0; Rscript Scripts/maftools.R {params.dir}/mutect2_out/oncotator_out/ mutect2_merged.maf {params.dir}/mutect2_out/mutect2_maf_summary {output.mutect2oncoprint}; Rscript Scripts/maftools.R {params.dir}/mutect_out/oncotator_out/ mutect_merged.maf {params.dir}/mutect_out/mutect_maf_summary {output.mutectoncoprint}; Rscript Scripts/maftools.R {params.dir}/strelka_out/oncotator_out/ strelka_merged.maf {params.dir}/strelka_out/strelka_maf_summary {output.strelkaoncoprint}"