rule vardict_tumoronly:
       input:  tumor="{x}.recal.bam",
               targets=config["project"]["workpath"]+"/exome_targets.bed",
               tumorbai="{x}.recal.bam.bai"
       output: vcf=config['project']['workpath']+"/vardict_out/{x}.vcf",
               filtvcf=config['project']['workpath']+"/vardict_out/{x}.FINAL.vcf",
               out=config['project']['workpath']+"/vardict_out/{x}.snpeff.out",
               csvstats=config['project']['workpath']+"/vardict_out/{x}.vardict.stats.csv",
               htmlstats=config['project']['workpath']+"/vardict_out/{x}.vardict.stats.html",
       params: sample="{x}",pon=config['references'][pfamily]['PON'],targets="exome_targets.bed",knowns=config['references'][pfamily]['MUTECTVARIANTS'],mutect=config['bin'][pfamily]['MUTECT'],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['MUTECTGENOME'],cosmic=config['references'][pfamily]['MUTECTCOSMIC'],snp=config['references'][pfamily]['MUTECTSNP'],snpeffgenome=config['references'][pfamily]['SNPEFF_GENOME'],snpeff=config['bin'][pfamily]['SNPEFF'],effconfig=config['references'][pfamily]['SNPEFF_CONFIG'],rname="pl:vardict"
       shell:  "mkdir -p vardict_out; module load R/3.5; module load samtools/1.9; /data/CCBR_Pipeliner/db/PipeDB/bin/VarDictJava/build/install/VarDict/bin/VarDict -G {params.genome} -f 0.05 -x 500 --nosv -b {input.tumor} -t -Q 20 -c 1 -S 2 -E 3 {params.targets} | /data/CCBR_Pipeliner/db/PipeDB/bin/VarDictJava/build/install/VarDict/bin/teststrandbias.R | /data/CCBR_Pipeliner/db/PipeDB/bin/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl -N {params.sample} -Q 20 -d 10 -v 6 -S -E -f 0.05 > {output.vcf}; module load GATK/3.8-0; module load java/1.8.0_92; GATK -m 48G SelectVariants -R {params.genome} --variant {output.vcf} --discordance {params.pon} -o {output.filtvcf}; module load snpEff/4.3t; java -Xmx12g -jar $SNPEFF_JAR -v {params.snpeffgenome} -c {params.effconfig} -cancer -canon -csvStats {output.csvstats} -stats {output.htmlstats} -cancerSamples pairs {output.filtvcf} > {output.out}"
