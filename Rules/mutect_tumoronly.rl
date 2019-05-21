rule mutect_tumoronly:
       input:  tumor="{x}.recal.bam",
               tumorbai="{x}.recal.bam.bai"
       output: vcf=config['project']['workpath']+"/mutect_out/{x}.vcf",
               stats=config['project']['workpath']+"/mutect_out/{x}.stats.out",
               vcfRename=config['project']['workpath']+"/mutect_out/{x}.FINAL.vcf",
               csvstats=config['project']['workpath']+"/mutect_out/{x}.mutect.stats.csv",
               htmlstats=config['project']['workpath']+"/mutect_out/{x}.mutect.stats.html",
               out=config['project']['workpath']+"/mutect_out/{x}.snpeff.out"
       params: knowns=config['references'][pfamily]['MUTECTVARIANTS'],mutect=config['bin'][pfamily]['MUTECT'],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['MUTECTGENOME'],cosmic=config['references'][pfamily]['MUTECTCOSMIC'],snp=config['references'][pfamily]['MUTECTSNP'],snpeffgenome=config['references'][pfamily]['SNPEFF_GENOME'],snpeff=config['bin'][pfamily]['SNPEFF'],effconfig=config['references'][pfamily]['SNPEFF_CONFIG'],rname="pl:mutect"
       shell:  "mkdir -p mutect_out; module load muTect/1.1.7; muTect --analysis_type MuTect --reference_sequence {params.genome} --vcf {output.vcf} {params.knowns} --disable_auto_index_creation_and_locking_when_reading_rods --input_file:tumor {input.tumor} --out {output.stats} -rf BadCigar; module load GATK/3.8-0; module load java/1.8.0_92; GATK -m 48G SelectVariants -R {params.genome} --variant {output.vcf} --excludeFiltered -o {output.vcfRename}; module load snpEff/4.3t; java -Xmx12g -jar $SNPEFF_JAR -v {params.snpeffgenome} -c {params.effconfig} -cancer -canon -csvStats {output.csvstats} -stats {output.htmlstats} -cancerSamples pairs {output.vcfRename} > {output.out}"