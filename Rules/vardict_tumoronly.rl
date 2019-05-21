rule vardict_tumoronly:
       input:  tumor="{x}.recal.bam",
               targets=ancient("exome_targets.bed"),
               tumorbai="{x}.recal.bam.bai"
       output: vcf=temp(config['project']['workpath']+"/vardict_out/{x}.vcf"),
               out=config['project']['workpath']+"/vardict_out/{x}.snpeff.out",
               csvstats=config['project']['workpath']+"/vardict_out/{x}.vardict.stats.csv",
               htmlstats=config['project']['workpath']+"/vardict_out/{x}.vardict.stats.html",
       params: targets="exome_targets.bed",knowns=config['references'][pfamily]['MUTECTVARIANTS'],mutect=config['bin'][pfamily]['MUTECT'],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['MUTECTGENOME'],cosmic=config['references'][pfamily]['MUTECTCOSMIC'],snp=config['references'][pfamily]['MUTECTSNP'],snpeffgenome=config['references'][pfamily]['SNPEFF_GENOME'],snpeff=config['bin'][pfamily]['SNPEFF'],effconfig=config['references'][pfamily]['SNPEFF_CONFIG'],rname="pl:vardict"
       shell:  "mkdir -p vardict_out; module load R/3.5; module load samtools/1.9; /data/CCBR_Pipeliner/db/PipeDB/bin/VarDict/vardict -G {params.genome} -f 0.05 -x 1000 --nosv -b {input.tumor} -y -t -Q 20 -c 1 -S 2 -E 3 {params.targets} | /data/CCBR_Pipeliner/db/PipeDB/bin/VarDict/teststrandbias.R | /data/CCBR_Pipeliner/db/PipeDB/bin/VarDict/var2vcf_valid.pl -S -E -f 0.05 2> {output.vcf}; module load snpEff/4.3t; java -Xmx12g -jar $SNPEFF_JAR -v {params.snpeffgenome} -c {params.effconfig} -cancer -canon -csvStats {output.csvstats} -stats {output.htmlstats} -cancerSamples pairs {output.vcf} > {output.out}"