rule alview:
      input: mutect="mutect_out/oncotator_out/mutect_variants.maf",
             strelka="strelka_out/oncotator_out/strelka_variants.maf"
      output: mutect="mutect_variants_alview.input",
              strelka="strelka_variants_alview.input"
      params: genome=config['references'][pfamily]['ALVIEWGENOME'],workdir=config['project']['workpath'],rname="pl:alview"
      shell:  "perl Scripts/run_alview.pl {input.mutect} {params.genome} mutect_out/oncotator_out/ mutect {params.workdir}; perl Scripts/run_alview.pl {input.strelka} {params.genome} strelka_out/oncotator_out/ strelka {params.workdir}; cd /data/CCBR_Pipeliner/db/PipeDB/bin/ALVIEW; ./alvgenslideshow mutect_variants_images {params.workdir}/{output.mutect} {params.genome}; ./alvgenslideshow strelka_variants_images {params.workdir}/{output.strelka} {params.genome}; mv mutect_variants_images.tar {params.workdir}/mutect_variants_images.tar; mv strelka_variants_images.tar {params.workdir}/strelka_variants_images.tar; rm -rf mutect_variants_images strelka_variants_images; rm mutect_variants_images.html strelka_variants_images.html"