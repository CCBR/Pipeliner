rule make_somatic_network:
    input: expand("germline_vcfs/{x}.vcf", x=samples)
    output: network="sample_network.bmp",
            vcf=temp("samples_and_knowns.vcf")
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['KNOWNANCESTRY'],regions="exome_targets.bed",rname="make.somatic.network"
    run:
      fl=os.popen("ls germline_vcfs/*.vcf").read().split()      
      var=" --variant "+" --variant ".join(fl)
      cmd="{params.gatk} -T CombineVariants -R {params.genome} --filteredrecordsmergetype KEEP_UNCONDITIONAL --genotypemergeoption UNIQUIFY -o {output.vcf} --variant {params.knowns} --minimumN 2"+var
      shell(cmd)
      cmd="module load vcftools; vcftools --vcf {output.vcf} --bed {params.regions} --recode --recode-INFO-all --out samples_and_knowns; perl Scripts/make_sample_network.pl samples_and_knowns.recode.vcf"
      shell(cmd)