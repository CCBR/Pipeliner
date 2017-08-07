rule vcf2maf_mutect:
  input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
         tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
         vcf=config['project']['workpath']+"/mutect_out/{x}.FINAL.vcf"
  output: maf=config['project']['workpath']+"/mutect_out/oncotator_out/{x}.maf"
  params: genome=config['references'][pfamily]['VEPGENOME'],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],build=config['references'][pfamily]['VEPBUILD'],species=config['references'][pfamily]['VEPSPECIES'],filtervcf=config['references'][pfamily]['VEPFILTERVCF'],rname="pl:vcf2maf"
  threads: 8
  shell: "module load vcf2maf/1.6.12; module load samtools/1.5; module load VEP/88; vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} --vep-path $VEP_HOME --vep-data $VEPCACHEDIR --ref-fasta {params.genome} --filter-vcf {params.filtervcf} --vep-forks 8 --vcf-tumor-id {params.tumorsample} --vcf-normal-id {params.normalsample} --tumor-id {params.tumorsample} --normal-id {params.normalsample} --ncbi-build {params.build} --species {params.species}"

rule vcf2maf_mutect2:
  input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
         tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
         vcf=config['project']['workpath']+"/mutect2_out/{x}.FINALmutect2.vcf"
  output: maf=config['project']['workpath']+"/mutect2_out/oncotator_out/{x}.maf"
  params: genome=config['references'][pfamily]['VEPGENOME'],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],build=config['references'][pfamily]['VEPBUILD'],species=config['references'][pfamily]['VEPSPECIES'],filtervcf=config['references'][pfamily]['VEPFILTERVCF'],rname="pl:vcf2maf"
  shell: "module load vcf2maf/1.6.12; module load samtools/1.5; module load VEP/88; vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} --vep-path $VEP_HOME --vep-data $VEPCACHEDIR --ref-fasta {params.genome} --filter-vcf {params.filtervcf} --vep-forks 8 --vcf-tumor-id {params.tumorsample} --vcf-normal-id {params.normalsample} --tumor-id {params.tumorsample} --normal-id {params.normalsample} --ncbi-build {params.build} --species {params.species}"

rule vcf2maf_strelka:
  input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
         tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
         vcf=config['project']['workpath']+"/strelka_out/{x}_filtered.vcf"
  output: maf=config['project']['workpath']+"/strelka_out/oncotator_out/{x}.maf"
  params: genome=config['references'][pfamily]['VEPGENOME'],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],build=config['references'][pfamily]['VEPBUILD'],species=config['references'][pfamily]['VEPSPECIES'],filtervcf=config['references'][pfamily]['VEPFILTERVCF'],rname="pl:vcf2maf"
  shell: "module load vcf2maf/1.6.12; module load samtools/1.5; module load VEP/88; vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} --vep-path $VEP_HOME --vep-data $VEPCACHEDIR --ref-fasta {params.genome} --filter-vcf {params.filtervcf} --vep-forks 8 --vcf-tumor-id {params.tumorsample} --vcf-normal-id {params.normalsample} --tumor-id {params.tumorsample} --normal-id {params.normalsample} --ncbi-build {params.build} --species {params.species}"