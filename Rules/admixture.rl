rule make_germline_network:
    input: "samples_and_knowns.vcf"
    output: dir="admixture",
            vcf="samples_and_knowns.vcf"
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['KNOWNANCESTRY'],rname="make.germline.network"
    shell: """
         {params.gatk} -T CombineVariants -R {params.genome} --genotypemergeoption UNSORTED -o {output.vcf} --variant {params.knowns} --variant {input} --minimumN 2; perl Scripts/make_sample_network.pl {output.vcf}

           """