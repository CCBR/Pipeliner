rule make_somatic_network:
    input: config['project']['workpath']+"/germline_vcfs/combined.vcf"
    output: network=config['project']['workpath']+"/sample_network_mqc.png",
    params: regions="exome_targets.bed",rname="make.somatic.network"
    run:
      cmd="module load vcftools; vcftools --vcf {input} --bed {params.regions} --recode --recode-INFO-all --out samples; module load perl/5.18.4; perl Scripts/make_sample_network.pl samples.recode.vcf; module load R; Rscript Scripts/magick.R; rm sample_network.bmp"
      shell(cmd)