rule make_genomeseq_network:
    input: "combined.strictFilter.vcf"
    output: network="sample_network.bmp",
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],rname="make.genomeseq.network"
    shell: """
         perl Scripts/make_sample_network.pl {input}

           """