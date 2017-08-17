if config['project']['annotation'] == "hg19":

  rule all-wgs-somatic-tumoronly:
    input:  expand("{s}"+".recal.bam",s=samples),
            expand("{s}"+".g.vcf",s=samples),
            expand(config['project']['workpath']+"/mutect2_out/{s}"+".FINALmutect2.vcf",s=samples),
            expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{s}"+".maf",s=samples),
            config['project']['workpath']+"/mutect2_out/mutect2_maf_summary.pdf",
            config['project']['workpath']+"/mutect2_out",
            config['project']['workpath']+"/cnvkit_out",
            config['project']['workpath']+"/delly_out",
            expand(config['project']['workpath']+"/cnvkit_out/{s}_calls.cns", s=samples),
            expand(config['project']['workpath']+"/cnvkit_out/{s}_gainloss.tsv", s=samples),                        
            config['project']['workpath']+"/delly_out/deletions.bcf",
            config['project']['workpath']+"/delly_out/inversions.bcf",
            config['project']['workpath']+"/delly_out/translocations.bcf",
            config['project']['workpath']+"/delly_out/duplications.bcf",
            config['project']['workpath']+"/delly_out/insertions.bcf",
            config['project']['workpath']+"/cnvkit_out/CNVkit_summary_heatmap.pdf",
            config['project']['workpath']+"/sample_network.bmp",
            config['project']['workpath']+"/mutect2_out/mutect2_variants.database",
            config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_variants.maf",
            config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
            expand("manta_out/{s}/results/variants/candidateSV.vcf.gz", s=samples),
            "admixture_out/admixture_table.tsv"
    output:
    params: rname="final"
    shell:  """
             module load multiqc/1.1; multiqc -f -e featureCounts .; mv *.out slurmfiles/; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped logfiles/

            """
elif config['project']['annotation'] == "hg38":

  rule all-wgs-somatic-tumoronly:
    input:  expand("{s}"+".recal.bam",s=samples),
            expand("{s}"+".g.vcf",s=samples),
            expand(config['project']['workpath']+"/mutect_out/{s}"+".FINAL.vcf",s=samples),
            expand(config['project']['workpath']+"/mutect2_out/{s}"+".FINALmutect2.vcf",s=samples),
            expand(config['project']['workpath']+"/strelka_out/{s}"+"_FINAL.vcf",s=samples),
            expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{s}"+".maf",s=samples),
            expand(config['project']['workpath']+"/strelka_out/oncotator_out/{s}"+".maf",s=samples),
            expand(config['project']['workpath']+"/mutect_out/oncotator_out/{s}"+".maf",s=samples),
            config['project']['workpath']+"/mutect2_out/mutect2_maf_summary.pdf",
            config['project']['workpath']+"/mutect_out/mutect_maf_summary.pdf",
            config['project']['workpath']+"/strelka_out/strelka_maf_summary.pdf",
#            config['project']['workpath']+"/mutect_variants_alview.input",
            config['project']['workpath']+"/strelka_out",
            config['project']['workpath']+"/mutect2_out",
            config['project']['workpath']+"/cnvkit_out",
            config['project']['workpath']+"/mutect_out",
            config['project']['workpath']+"/delly_out",
            expand(config['project']['workpath']+"/cnvkit_out/{s}_calls.cns", s=samples),
            expand(config['project']['workpath']+"/cnvkit_out/{s}_gainloss.tsv", s=samples),                        
            expand(config['project']['workpath']+"/delly_out/{s}_del.bcf", s=samples),
            expand(config['project']['workpath']+"/delly_out/{s}_ins.bcf", s=samples),
            expand(config['project']['workpath']+"/delly_out/{s}_dup.bcf", s=samples),
            expand(config['project']['workpath']+"/delly_out/{s}_tra.bcf", s=samples),
            expand(config['project']['workpath']+"/delly_out/{s}_inv.bcf", s=samples),
            expand(config['project']['workpath']+"/theta_out/{s}/{s}_thetaIN", s=samples),
            expand(config['project']['workpath']+"/conpair_out/{s}.conpair", s=samples),
            config['project']['workpath']+"/cnvkit_out/CNVkit_summary_heatmap.pdf",
            config['project']['workpath']+"/sample_network.bmp",
            config['project']['workpath']+"/mutect2_out/mutect2_variants.database",
            config['project']['workpath']+"/mutect_out/mutect_variants.database",
            config['project']['workpath']+"/strelka_out/strelka_variants.database",
            config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_variants.maf",
            config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/strelka_out/oncotator_out/strelka_variants.maf",
            config['project']['workpath']+"/strelka_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/mutect_out/oncotator_out/mutect_variants.maf",
            config['project']['workpath']+"/mutect_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/exome_targets.bed",
            expand("manta_out/{s}/results/variants/candidateSV.vcf.gz", s=samples),
            "admixture_out/admixture_table.tsv"
    output:
    params: rname="final"
    shell:  """
             module load multiqc/1.1; multiqc -f -e featureCounts .; mv *.out slurmfiles/; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped logfiles/

            """

elif config['project']['annotation'] == "mm10":

  rule all-wgs-somatic-tumoronly:
    input:  expand("{s}"+".recal.bam",s=samples),
            expand(config['project']['workpath']+"/mutect_out/{s}"+".FINAL.vcf",s=samples),
            expand(config['project']['workpath']+"/mutect2_out/{s}"+".FINALmutect2.vcf",s=samples),
            expand(config['project']['workpath']+"/strelka_out/{s}"+"_FINAL.vcf",s=samples),
            expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{s}"+".maf",s=samples),
            expand(config['project']['workpath']+"/strelka_out/oncotator_out/{s}"+".maf",s=samples),
            expand(config['project']['workpath']+"/mutect_out/oncotator_out/{s}"+".maf",s=samples),
#            config['project']['workpath']+"/mutect_variants_alview.input",
            config['project']['workpath']+"/strelka_out",
            config['project']['workpath']+"/mutect2_out",
            config['project']['workpath']+"/cnvkit_out",
            config['project']['workpath']+"/mutect_out",
            config['project']['workpath']+"/delly_out",
            expand(config['project']['workpath']+"/cnvkit_out/{s}_calls.cns", s=samples),
            expand(config['project']['workpath']+"/cnvkit_out/{s}_gainloss.tsv", s=samples),                        
            expand(config['project']['workpath']+"/delly_out/{s}_del.bcf", s=samples),
            expand(config['project']['workpath']+"/delly_out/{s}_ins.bcf", s=samples),
            expand(config['project']['workpath']+"/delly_out/{s}_dup.bcf", s=samples),
            expand(config['project']['workpath']+"/delly_out/{s}_tra.bcf", s=samples),
            expand(config['project']['workpath']+"/delly_out/{s}_inv.bcf", s=samples),
            expand(config['project']['workpath']+"/theta_out/{s}/{s}_thetaIN", s=samples),
            expand(config['project']['workpath']+"/conpair_out/{s}.conpair", s=samples),
            config['project']['workpath']+"/cnvkit_out/CNVkit_summary_heatmap.pdf",
            config['project']['workpath']+"/mutect2_out/mutect2_variants.database",
            config['project']['workpath']+"/mutect_out/mutect_variants.database",
            config['project']['workpath']+"/strelka_out/strelka_variants.database",
            config['project']['workpath']+"/sample_network.bmp",
            config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_variants.maf",
            config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/strelka_out/oncotator_out/strelka_variants.maf",
            config['project']['workpath']+"/strelka_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/mutect_out/oncotator_out/mutect_variants.maf",
            config['project']['workpath']+"/mutect_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/exome_targets.bed",
            expand("manta_out/{s}/results/variants/candidateSV.vcf.gz", s=samples),
            "admixture_out/admixture_table.tsv"
    output:
    params: rname="final"
    shell:  """
             module load multiqc/1.1; multiqc -f -e featureCounts .; mv *.out slurmfiles/; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped logfiles/

            """