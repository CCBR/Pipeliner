rule gatk_apply_vqsr:
	input: vcf='../gvcf/all.vcf', recal='../gvcf/all.{mode}.recal', tranches ='../gvcf/all.{mode}.tranches'
	output: vcff='../gvcf/all.vqsr.{mode}.vcf'
	params: batch="-l nodes=1:gpfs", GATK_BASE=config.GATK_BASE2, REF=config.REF, Java=config.JAVA,,rname="pl:vqsr"
	run: 
		shell ('{params.Java}/{params.GATK_BASE} -T ApplyRecalibration -R {params.REF} -input {input.vcf} --ts_filter_level 99.0 -tranchesFile {input.tranches} -recalFile {input.recal} -mode {wildcards.mode} -o {output.vcff}')
