rule alview:
      input: "mutect_out/oncotator_out/mutect_variants.maf"
      output: "mutect_out/oncotator_out/mutect_variants_alview.input"
      params: genome=config['references'][pfamily]['ALVIEWGENOME'],workdir=config['project']['workpath'],rname="pl:alview"
      shell:  "perl Scripts/run_alview.pl {input} {params.genome} {params.workdir}; cd /data/CCBR/apps/ALVIEW; ./alvgenslideshow mutect_variants_images {params.workdir}/mutect_variants_alview.input {params.genome}; mv mutect_variants_images.tar {params.workdir}/mutect_variants_images.tar; rm -rf mutect_variants_images; rm mutect_variants_images.html"