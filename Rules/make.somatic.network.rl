rule make_somatic_network:
    input: config['project']['workpath']+"/germline_vcfs/combined.vcf"
    output: network=config['project']['workpath']+"/sample_network.bmp",
            vcf=temp(config['project']['workpath']+"/samples_and_knowns.vcf")
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['KNOWNANCESTRY'],regions="exome_targets.bed",rname="make.somatic.network"
    run:
      cmd="{params.gatk} -T CombineVariants -R {params.genome} --filteredrecordsmergetype KEEP_UNCONDITIONAL --genotypemergeoption UNSORTED -o {output.vcf} --variant {params.knowns} --variant {input} --minimumN 2"
      shell(cmd)
      cmd="module load vcftools; vcftools --vcf {output.vcf} --bed {params.regions} --recode --recode-INFO-all --out samples_and_knowns; perl Scripts/make_sample_network.pl samples_and_knowns.recode.vcf"
      shell(cmd)