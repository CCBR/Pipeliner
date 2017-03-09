rule alview:
      input: mutect=config['project']['workpath']+"/mutect_out/oncotator_out/mutect_variants.maf",
             strelka=config['project']['workpath']+"/strelka_out/oncotator_out/strelka_variants.maf",
             mutect2=config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_variants.maf"
      output: mutect=config['project']['workpath']+"/mutect_variants_alview.input",
              strelka=config['project']['workpath']+"/strelka_variants_alview.input",
              mutect2=config['project']['workpath']+"/mutect2_variants_alview.input"
      params: genome=config['references'][pfamily]['ALVIEWGENOME'],workdir=config['project']['workpath'],rname="pl:alview"
      shell:  "perl Scripts/run_alview.pl {input.mutect} {params.genome} mutect_out/oncotator_out/ mutect {params.workdir}; perl Scripts/run_alview.pl {input.strelka} {params.genome} strelka_out/oncotator_out/ strelka {params.workdir}; perl Scripts/run_alview.pl {input.mutect2} {params.genome} mutect2_out/oncotator_out/ mutect2 {params.workdir}; cd /data/CCBR_Pipeliner/db/PipeDB/bin/ALVIEW; ./alvgenslideshow {params.workdir}/mutect_variants_images {output.mutect} {params.genome}; ./alvgenslideshow {params.workdir}/strelka_variants_images {output.strelka} {params.genome}; ./alvgenslideshow {params.workdir}/mutect2_variants_images {output.mutect2} {params.genome}"