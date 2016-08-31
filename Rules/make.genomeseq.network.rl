rule make_genomeseq_network:
    input: "combined.vcf"
    output: network="sample_network.bmp",
            vcf="samples_and_knowns.vcf"
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['KNOWNANCESTRY'],rname="make.genomeseq.network"
    shell: """
         {params.gatk} -T CombineVariants -R {params.genome} --genotypemergeoption UNSORTED -o {output.vcf} --variant {params.knowns} --variant {input} --minimumN 2; perl Scripts/make_sample_network.pl {output.vcf}

           """