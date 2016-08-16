rule alview:
      input: "mutect_out/oncotator_out/mutect_variants.maf"
      output: "mutect_out/oncotator_out/mutect_variants_alview.input"
      params: genome=config['references'][pfamily]['ALVIEWGENOME'],rname="pl:alview"
      shell:  "perl Scripts/run_alview.pl {input} {params.genome}"