rule make_germline_network:
    input: "exome.recode.vcf"
    output:"sample_network.bmp"
    params: batch ="-l nodes=1:gpfs -q ccr",rname="make.germline.network"
    shell: """
         perl Scripts/make_sample_network.pl {input}

           """