rule make_somatic_network:
    input: vcf="germline_snps.vcf"
    output: network="sample_network.bmp",
            vcf="samples_and_knowns.vcf"
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],rname="make.somatic.network"
    shell: """
         perl Scripts/make_sample_network.pl {output.vcf}

           """