rule mutect:
       input:  normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
               tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
               targets=ancient("exome_targets.bed")
       output: vcf=config['project']['workpath']+"/mutect_out/{x}.vcf",
               stats=config['project']['workpath']+"/mutect_out/{x}.stats.out",
               vcfRename=config['project']['workpath']+"/mutect_out/{x}.FINAL.vcf",
               csvstats=config['project']['workpath']+"/mutect_out/{x}.mutect.stats.csv",
               htmlstats=config['project']['workpath']+"/mutect_out/{x}.mutect.stats.html",
               out=config['project']['workpath']+"/mutect_out/{x}.snpeff.out"
       params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],targets="exome_targets.bed",knowns=config['references'][pfamily]['MUTECTVARIANTS'],mutect=config['bin'][pfamily]['MUTECT'],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['MUTECTGENOME'],cosmic=config['references'][pfamily]['MUTECTCOSMIC'],snp=config['references'][pfamily]['MUTECTSNP'],snpeffgenome=config['references'][pfamily]['SNPEFF_GENOME'],snpeff=config['bin'][pfamily]['SNPEFF'],rname="pl:mutect"
       shell:  "{params.mutect} --analysis_type MuTect --reference_sequence {params.genome} --vcf {output.vcf} {params.knowns} --intervals {params.targets} --input_file:normal {input.normal} --input_file:tumor {input.tumor} --out {output.stats} -rf BadCigar; {params.gatk} -T SelectVariants -R {params.genome} --variant {output.vcf} --excludeFiltered -o {output.vcfRename}; module load snpEff; {params.snpeff} -v {params.snpeffgenome} -interval {params.targets} -cancer -canon -csvStats {output.csvstats} -stats {output.htmlstats} -cancerSamples pairs {output.vcfRename} > {output.out}"