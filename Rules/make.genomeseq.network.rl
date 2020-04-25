rule make_genomeseq_network:
    input: "combined.strictFilter.vcf"
    output: network="sample_network_mqc.png",
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],rname="make.genomeseq.network"
    shell: """
         module load perl/5.18.4; perl Scripts/make_sample_network.pl {input}; module load R; Rscript Scripts/magick.R; rm sample_network.bmp

           """