rule make_wgs_somatic_network:
    input: config['project']['workpath']+"/germline_vcfs/combined.vcf"
    output: network=config['project']['workpath']+"/sample_network_mqc.png",
    params: rname="make.somatic.network"
    run:
      cmd="module load vcftools; vcftools --vcf {input} --recode --recode-INFO-all --out samples; module load perl/5.18.4; perl Scripts/make_sample_network.pl samples.recode.vcf; module load R; Rscript Scripts/magick.R; rm sample_network.bmp"
      shell(cmd)