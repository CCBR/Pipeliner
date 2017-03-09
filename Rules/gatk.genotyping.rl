rule gatk_genotyping:
	input: gvcf = expand("../gvcf/{sample}.gvcf", sample=samples)
	output: '../gvcf/all.vcf'
	params: batch="-l nodes=1:gpfs", GATK_BASE=config.GATK_BASE2, REF=config.REF, snp138=config.SNP138,rname="pl:genotyping"
        threads: 4
        run: gvcf = gatk_multi_arg(" --variant", input.gvcf);shell('/usr/local/java/jdk1.7.0_02/bin/{params.GATK_BASE} -T GenotypeGVCFs -R {params.REF} --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts {gvcf} --dbsnp {params.snp138} -o {output}')
