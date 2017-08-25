rule make_germline_network:
    input: "exome.strictFilter.vcf"
    output: network="sample_network.bmp",
    params: rname="make.germline.network"
    shell: """
         perl Scripts/make_sample_network.pl {input}

           """