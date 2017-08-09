rule make_genomeseq_network:
    input: "combined.strictFilter.vcf"
    output: network="sample_network.bmp",
            vcf="samples_and_knowns.vcf"
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['KNOWNANCESTRY'],rname="make.genomeseq.network"
    shell: """
         perl Scripts/make_sample_network.pl {output.vcf}

           """