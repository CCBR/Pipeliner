if config['project']['annotation'] == "hg19":

  rule all_exomeseq_somatic:
    input:  expand("{s}"+".recal.bam",s=samples),
            expand("{s}"+".g.vcf",s=samples),
            expand(config['project']['workpath']+"/mutect_out/{p}"+".FINAL.vcf",p=pairs),
            expand(config['project']['workpath']+"/mutect2_out/{p}"+".FINALmutect2.vcf",p=pairs),
            expand(config['project']['workpath']+"/strelka_out/{p}"+".vcf",p=pairs),
            expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/strelka_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/mutect_out/oncotator_out/{p}"+".maf",p=pairs),
#            config['project']['workpath']+"/mutect_variants_alview.input",
            config['project']['workpath']+"/strelka_out",
            config['project']['workpath']+"/mutect2_out",
            config['project']['workpath']+"/cnvkit_out",
            config['project']['workpath']+"/mutect_out",
            config['project']['workpath']+"/delly_out",
            expand(config['project']['workpath']+"/cnvkit_out/{p}_calls.cns", p=pairs),
            expand(config['project']['workpath']+"/cnvkit_out/{p}_gainloss.tsv", p=pairs),                        
            expand(config['project']['workpath']+"/delly_out/{p}_del.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_ins.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_dup.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_tra.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_inv.bcf", p=pairs),
            expand(config['project']['workpath']+"/theta_out/{p}/{p}_thetaIN", p=pairs),
            expand(config['project']['workpath']+"/conpair_out/{p}.conpair", p=pairs),
            config['project']['workpath']+"/mutect_out/merged_somatic.vcf",
            config['project']['workpath']+"/cnvkit_out/CNVkit_summary_heatmap.pdf",
            config['project']['workpath']+"/mutect_out/merged_somatic_snpEff.vcf",
            config['project']['workpath']+"/strelka_out/merged_somatic_snpEff.vcf",
#            config['project']['workpath']+"/mutect2_out/merged_somatic.vcf",
#            config['project']['workpath']+"/mutect2_out/merged_somatic_snpEff.vcf",
#            config['project']['workpath']+"/variants.bed",
#            config['project']['workpath']+"/full_annot.txt.zip",
            config['project']['workpath']+"/sample_network.bmp",
#            config['project']['workpath']+"/variants.database",
            config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_variants.maf",
            config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/strelka_out/oncotator_out/strelka_variants.maf",
            config['project']['workpath']+"/strelka_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/mutect_out/oncotator_out/mutect_variants.maf",
            config['project']['workpath']+"/mutect_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/exome_targets.bed",
#            expand(config['project']['workpath']+"/manta_out/{p}/candidateSV.vcf.gz", p=pairs)
    output:
    params: rname="final"
    shell:  """
             module load multiqc; multiqc -f -e featureCounts .; mv *.out slurmfiles/; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped logfiles/

            """

else:

  rule all_exomeseq_somatic:
    input:  expand("{s}"+".recal.bam",s=samples),
            expand(config['project']['workpath']+"/mutect_out/{p}"+".FINAL.vcf",p=pairs),
            expand(config['project']['workpath']+"/mutect2_out/{p}"+".FINALmutect2.vcf",p=pairs),
            expand(config['project']['workpath']+"/strelka_out/{p}"+".vcf",p=pairs),
            expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/strelka_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/mutect_out/oncotator_out/{p}"+".maf",p=pairs),
#            config['project']['workpath']+"/mutect_variants_alview.input",
            config['project']['workpath']+"/strelka_out",
            config['project']['workpath']+"/mutect2_out",
            config['project']['workpath']+"/cnvkit_out",
            config['project']['workpath']+"/mutect_out",
            config['project']['workpath']+"/delly_out",
            expand(config['project']['workpath']+"/cnvkit_out/{p}_calls.cns", p=pairs),
            expand(config['project']['workpath']+"/cnvkit_out/{p}_gainloss.tsv", p=pairs),                        
            expand(config['project']['workpath']+"/delly_out/{p}_del.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_ins.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_dup.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_tra.bcf", p=pairs),
            expand(config['project']['workpath']+"/delly_out/{p}_inv.bcf", p=pairs),
            expand(config['project']['workpath']+"/theta_out/{p}/{p}_thetaIN", p=pairs),
            expand(config['project']['workpath']+"/conpair_out/{p}.conpair", p=pairs),
            config['project']['workpath']+"/mutect_out/merged_somatic.vcf",
            config['project']['workpath']+"/cnvkit_out/CNVkit_summary_heatmap.pdf",
            config['project']['workpath']+"/mutect_out/merged_somatic_snpEff.vcf",
            config['project']['workpath']+"/strelka_out/merged_somatic_snpEff.vcf",
#            config['project']['workpath']+"/mutect2_out/merged_somatic.vcf",
#            config['project']['workpath']+"/mutect2_out/merged_somatic_snpEff.vcf",
#            config['project']['workpath']+"/variants.bed",
#            config['project']['workpath']+"/full_annot.txt.zip",
            config['project']['workpath']+"/sample_network.bmp",
#            config['project']['workpath']+"/variants.database",
            config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_variants.maf",
            config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/strelka_out/oncotator_out/strelka_variants.maf",
            config['project']['workpath']+"/strelka_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/mutect_out/oncotator_out/mutect_variants.maf",
            config['project']['workpath']+"/mutect_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/exome_targets.bed",
#            expand(config['project']['workpath']+"/manta_out/{p}/candidateSV.vcf.gz", p=pairs)
    output:
    params: rname="final"
    shell:  """
             module load multiqc; multiqc -f -e featureCounts .; mv *.out slurmfiles/; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped logfiles/

            """