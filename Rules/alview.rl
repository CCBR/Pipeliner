rule alview:
      input: mutect=config['project']['workpath']+"/mutect_out/oncotator_out/mutect_variants.maf",
             strelka=config['project']['workpath']+"/strelka_out/oncotator_out/strelka_variants.maf"
      output: mutect=config['project']['workpath']+"/mutect_variants_alview.input",
              strelka=config['project']['workpath']+"/strelka_variants_alview.input"
      params: genome=config['references'][pfamily]['ALVIEWGENOME'],workdir=config['project']['workpath'],rname="pl:alview"
      shell:  "perl Scripts/run_alview.pl {input.mutect} {params.genome} mutect_out/oncotator_out/ mutect {params.workdir}; perl Scripts/run_alview.pl {input.strelka} {params.genome} strelka_out/oncotator_out/ strelka {params.workdir}; cd /data/CCBR_Pipeliner/db/PipeDB/bin/ALVIEW; ./alvgenslideshow {params.workdir}/mutect_variants_images {params.workdir}/{output.mutect} {params.genome}; ./alvgenslideshow {params.workdir}/strelka_variants_images {params.workdir}/{output.strelka} {params.genome}"