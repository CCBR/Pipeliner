rule maftools:
  input: expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{x}.maf",x=pairs),
         expand(config['project']['workpath']+"/strelka_out/oncotator_out/{x}.maf",x=pairs),
         expand(config['project']['workpath']+"/mutect_out/oncotator_out/{x}.maf",x=pairs),
  output: mutectpre=config['project']['workpath']+"/mutect_out/oncotator_out/mutect_variants.maf",
          mutect2pre=config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_variants.maf",
          strelkapre=config['project']['workpath']+"/strelka_out/oncotator_out/strelka_variants.maf",
          mutect2maf=config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_filtered.maf",
          mutectmaf=config['project']['workpath']+"/mutect_out/oncotator_out/mutect_filtered.maf",
          strelkamaf=config['project']['workpath']+"/strelka_out/oncotator_out/strelka_filtered.maf",
          mutect2summary=config['project']['workpath']+"/mutect2_out/mutect2_maf_summary.pdf",
          mutect2oncoprint=config['project']['workpath']+"/mutect2_out/mutect2_oncoplot.pdf",
          mutectsummary=config['project']['workpath']+"/mutect_out/mutect_maf_summary.pdf",
          mutectoncoprint=config['project']['workpath']+"/mutect_out/mutect_oncoplot.pdf",
          strelkasummary=config['project']['workpath']+"/strelka_out/strelka_maf_summary.pdf",
          strelkaoncoprint=config['project']['workpath']+"/strelka_out/strelka_oncoplot.pdf",
  params: dir=config['project']['workpath'],rname="pl:maftools"
  shell: "cat mutect2_out/oncotator_out/*.maf > mutect2_out/oncotator_out/mutect2_variants.maf; cat mutect_out/oncotator_out/*.maf > mutect_out/oncotator_out/mutect_variants.maf; cat strelka_out/oncotator_out/*.maf > strelka_out/oncotator_out/strelka_variants.maf; perl Scripts/prep_mafs.pl mutect2_out/oncotator_out/mutect2_variants.maf mutect2; perl Scripts/prep_mafs.pl mutect_out/oncotator_out/mutect_variants.maf mutect; perl Scripts/prep_mafs.pl strelka_out/oncotator_out/strelka_variants.maf strelka; module load R/3.4.0_gcc-6.2.0; Rscript Scripts/maftools.R {params.dir}/mutect2_out/oncotator_out/ mutect2_filtered.maf {params.dir}/mutect2_out/mutect2_maf_summary {output.mutect2oncoprint}; Rscript Scripts/maftools.R {params.dir}/mutect_out/oncotator_out/ mutect_filtered.maf {params.dir}/mutect_out/mutect_maf_summary {output.mutectoncoprint}; Rscript Scripts/maftools.R {params.dir}/strelka_out/oncotator_out/ strelka_filtered.maf {params.dir}/strelka_out/strelka_maf_summary {output.strelkaoncoprint}"