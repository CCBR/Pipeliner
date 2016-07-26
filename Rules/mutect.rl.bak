rule mutect:
       input:  normal=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".realign.bam",
               tumor=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".realign.bam"
       output: vcf="mutect/{x}.vcf",
               stats="mutect/{x}.stats.out",
               wig="mutect/{x}.wig.txt",
               vcfFilter="mutect/{x}.PASS.vcf",
               vcfRename="mutect/{x}.FINAL.vcf",
               names="mutect/{x}.samples"               
       params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],mutect=config['bin']['MUTECT'],gatk=config['bin']['GATK'],genome=config['references']['MUTECTGENOME'],cosmic=config['references']['MUTECTCOSMIC'],snp=config['references']['MUTECTSNP'],rname="pl:mutect"
       shell:  """
{params.mutect} --analysis_type MuTect --reference_sequence {params.genome} --vcf {output.vcf} --cosmic {params.cosmic} {params.snp} --input_file:normal {input.normal} --input_file:tumor {input.tumor} --out {output.stats} -rf BadCigar; {params.gatk} -T SelectVariants -R {params.genome} --variant {output.vcf} --excludeFiltered -o {output.vcfFilter}; gzip {output.vcfFilter}; echo $'{params.normalsample}\n{params.tumorsample}' > {output.names} ; module load samtools; bcftools reheader -o {output.vcfRename} -s {output.names}  {output.vcfFilter}.gz; gunzip {output.vcfFilter}.gz 

"""
#--coverage_file {output.wig}
#--dbsnp {params.snp}
