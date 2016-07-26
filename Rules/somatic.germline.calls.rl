rule somatic_germline_calls:
     input:  expand("{s}"+".realign.bam",s=samples)
     output: unfiltered="germline.vcf",
             filtered="germline_snps.vcf"
     params: genome=config['references'][pfamily]['GENOME'],regions=config['references'][pfamily]['REFFLAT'],rname="pl:germcalls"
     shell:  "ls *realign.bam > samtoolsinbams; module load samtools; samtools mpileup -ug -b samtoolsinbams --fasta-ref {params.genome} --positions {params.regions} | bcftools call -mvO v - > {output.unfiltered}; module load vcftools; vcftools --vcf {output.unfiltered} --minQ 20 --remove-indels --recode --out germsnps; mv germsnps.recode.vcf {output.filtered}"