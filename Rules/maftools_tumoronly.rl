rule maftools_tumoronly:
  input: mutect2=config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_variants.maf",
  output: mutect2maf=config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_merged.maf",
          mutect2summary=config['project']['workpath']+"/mutect2_out/mutect2_maf_summary.pdf",
          mutect2oncoprint=config['project']['workpath']+"/mutect2_out/mutect2_oncoplot.pdf",
  params: dir=config['project']['workpath'],rname="pl:maftools"
  shell: "perl Scripts/prep_mafs.pl {input.mutect2} mutect2; module load R/3.4.0_gcc-6.2.0; Rscript Scripts/maftools.R {params.dir}/mutect2_out/oncotator_out/ mutect2_merged.maf {params.dir}/mutect2_out/mutect2_maf_summary {output.mutect2oncoprint}"