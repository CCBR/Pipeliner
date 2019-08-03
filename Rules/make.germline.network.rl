rule make_germline_network:
    input: "exome.strictFilter.vcf.gz"
    output: network="sample_network_mqc.png",
    params: rname="make.germline.network"
    shell: """
         perl Scripts/make_sample_network.pl {input}; module load R; Rscript Scripts/magick.R; rm sample_network.bmp

           """