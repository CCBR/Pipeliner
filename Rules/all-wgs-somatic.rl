if config['project']['annotation'] == "hg19":
  rule all_wgs_somatic:
    input: "combined.vcf",
           config['project']['workpath']+"/full_annot.txt.zip",
#           "variants.database",
           "sample_network.bmp",
           "combined.snpeff.vcf",
           "combined.strictFilter.vcf",
            expand(config['project']['workpath']+"/delly_out/{p}_del.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_ins.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_dup.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_tra.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_inv.bcf", p=pairs),
            expand(config['project']['workpath']+"/theta_out/{p}/{p}_thetaIN", p=pairs),
            expand(config['project']['workpath']+"/conpair_out/{p}.conpair", p=pairs),
            expand("manta_out/{p}/results/variants/candidateSV.vcf.gz", p=pairs),
            expand(config['project']['workpath']+"/mutect_out/{p}"+".FINAL.vcf",p=pairs),
            expand(config['project']['workpath']+"/mutect2_out/{p}"+".FINALmutect2.vcf",p=pairs),
            expand(config['project']['workpath']+"/strelka_out/{p}"+".vcf",p=pairs),
            expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/strelka_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/mutect_out/oncotator_out/{p}"+".maf",p=pairs),
            config['project']['workpath']+"/strelka_out",
            config['project']['workpath']+"/mutect2_out",
            config['project']['workpath']+"/cnvkit_out",
            config['project']['workpath']+"/mutect_out",
            config['project']['workpath']+"/delly_out",
           "cnvkit_out/cnvkit_heatmap.pdf",
#           "breakdancer_out/file.ctx",
           config['project']['workpath']+"/svaba_out/{p}.log",
           "admixture_out/samples_and_knowns_filtered_recode.Q"
    output: 
    params: rname="final"
    shell:  """
             module load multiqc/1.1; multiqc -f -e featureCounts .; mv *.out slurmfiles/; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """

else:

  rule all_wgs_somatic:
    input: "combined.vcf",
           config['project']['workpath']+"/full_annot.txt.zip",
#           "variants.database",
           "sample_network.bmp",
           expand("sample_vcfs/{s}"+".stats.csv",s=samples),           
           "combined.snpeff.vcf",
           "combined.strictFilter.vcf",
           "manta_out/results/variants/diploidSV.vcf.gz",
            expand(config['project']['workpath']+"/delly_out/{p}_del.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_ins.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_dup.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_tra.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_inv.bcf", p=pairs),
            expand(config['project']['workpath']+"/theta_out/{p}/{p}_thetaIN", p=pairs),
            expand(config['project']['workpath']+"/conpair_out/{p}.conpair", p=pairs),
            expand("manta_out/{p}/results/variants/candidateSV.vcf.gz", p=pairs),
            expand(config['project']['workpath']+"/mutect_out/{p}"+".FINAL.vcf",p=pairs),
            expand(config['project']['workpath']+"/mutect2_out/{p}"+".FINALmutect2.vcf",p=pairs),
            expand(config['project']['workpath']+"/strelka_out/{p}"+".vcf",p=pairs),
            expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/strelka_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/mutect_out/oncotator_out/{p}"+".maf",p=pairs),
            config['project']['workpath']+"/strelka_out",
            config['project']['workpath']+"/mutect2_out",
            config['project']['workpath']+"/cnvkit_out",
            config['project']['workpath']+"/mutect_out",
            config['project']['workpath']+"/delly_out",
           "cnvkit_out/cnvkit_heatmap.pdf",
#           "breakdancer_out/file.ctx",
           config['project']['workpath']+"/svaba_out/{p}.log",
           "admixture_out/samples_and_knowns_filtered_recode.Q"
    output: 
    params: rname="final"
    shell:  """
             module load multiqc/1.1; multiqc -f -e featureCounts .; mv *.out slurmfiles/; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """