rule maftools_mutect2_tumoronly:
  input: mutect2mafs=expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{s}.maf",s=samples),
  output: pre=temp(config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_variants.maf"),
          maf=config['project']['workpath']+"/mutect2_out/oncotator_out/final_filtered.maf",
          summary=config['project']['workpath']+"/mutect2_out/mutect2_maf_summary.pdf",
          oncoprint=config['project']['workpath']+"/mutect2_out/mutect2_oncoplot.pdf",
  params: dir=config['project']['workpath'],rname="pl:maftools"
  shell: "cat {input} > mutect2_out/oncotator_out/mutect2_variants.maf; module load perl/5.18.4; perl Scripts/prep_mafs.pl mutect2_out/oncotator_out/mutect2_variants.maf mutect2_out; module load R/3.5; Rscript Scripts/maftools.R {params.dir}/mutect2_out/oncotator_out/ final_filtered.maf {params.dir}/mutect2_out/mutect2_maf_summary.pdf {output.oncoprint}"

rule maftools_tumoronly_mutect:
  input: mutectmaf=expand(config['project']['workpath']+"/mutect_out/oncotator_out/{s}.maf",s=samples),
  output: pre=temp(config['project']['workpath']+"/mutect_out/oncotator_out/mutect_variants.maf"),
          maf=config['project']['workpath']+"/mutect_out/oncotator_out/final_filtered.maf",
          summary=config['project']['workpath']+"/mutect_out/mutect_maf_summary.pdf",
          oncoprint=config['project']['workpath']+"/mutect_out/mutect_oncoplot.pdf",
  params: dir=config['project']['workpath'],rname="pl:maftools"
  shell: "cat {input} > mutect_out/oncotator_out/mutect_variants.maf; perl Scripts/prep_mafs.pl mutect_out/oncotator_out/mutect_variants.maf mutect_out; module load R/3.5; Rscript Scripts/maftools.R {params.dir}/mutect_out/oncotator_out/ final_filtered.maf {params.dir}/mutect_out/mutect_maf_summary.pdf {output.oncoprint}"

rule maftools_tumoronly_vardict:
  input: vardictmaf=expand(config['project']['workpath']+"/vardict_out/oncotator_out/{s}.maf",s=samples),
  output: pre=temp(config['project']['workpath']+"/vardict_out/oncotator_out/vardict_variants.maf"),
          maf=config['project']['workpath']+"/vardict_out/oncotator_out/final_filtered.maf",
          summary=config['project']['workpath']+"/vardict_out/vardict_maf_summary.pdf",
          oncoprint=config['project']['workpath']+"/vardict_out/vardict_oncoplot.pdf",
  params: dir=config['project']['workpath'],rname="pl:maftools"
  shell: "cat {input} > vardict_out/oncotator_out/vardict_variants.maf; perl Scripts/prep_mafs.pl vardict_out/oncotator_out/vardict_variants.maf vardict_out; module load R/3.5; Rscript Scripts/maftools.R {params.dir}/vardict_out/oncotator_out/ final_filtered.maf {params.dir}/vardict_out/vardict_maf_summary.pdf {output.oncoprint}"

rule maftools_tumoronly_merged:
  input: mergedmaf=expand(config['project']['workpath']+"/merged_somatic_variants/oncotator_out/{s}.maf",s=samples),
  output: pre=temp(config['project']['workpath']+"/merged_somatic_variants/oncotator_out/merged_variants.maf"),
          maf=config['project']['workpath']+"/merged_somatic_variants/oncotator_out/final_filtered.maf",
          summary=config['project']['workpath']+"/merged_somatic_variants/merged_maf_summary.pdf",
          oncoprint=config['project']['workpath']+"/merged_somatic_variants/merged_oncoplot.pdf",
  params: dir=config['project']['workpath'],rname="pl:maftools"
  shell: "cat {input} > merged_somatic_variants/oncotator_out/merged_variants.maf; perl Scripts/prep_mafs.pl merged_somatic_variants/oncotator_out/merged_variants.maf merged_somatic_variants; module load R/3.5; Rscript Scripts/maftools.R {params.dir}/merged_somatic_variants/oncotator_out/ final_filtered.maf {params.dir}/merged_somatic_variants/merged_maf_summary.pdf {output.oncoprint}"