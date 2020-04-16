if config['project']['annotation'] == "hg19":

  rule all_exomeseq_somatic_tumoronly:
    input:  expand("{s}"+".recal.bam",s=samples),
            expand(config['project']['workpath']+"/mutect2_out/{s}"+".FINALmutect2.vcf",s=samples),
            expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{s}"+".maf",s=samples),
            config['project']['workpath']+"/mutect2_out/mutect2_maf_summary.pdf",
            config['project']['workpath']+"/mutect_out/mutect_maf_summary.pdf",
            config['project']['workpath']+"/vardict_out/vardict_maf_summary.pdf",
            config['project']['workpath']+"/merged_somatic_variants/merged_maf_summary.pdf",
            config['project']['workpath']+"/mutect2_out",
            config['project']['workpath']+"/mutect2_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/mutect_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/vardict_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/merged_somatic_variants/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/mutect_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/vardict_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/merged_somatic_variants/mutsigCV_out/somatic.sig_genes.txt",
            expand(config['project']['workpath']+"/mutect_out/{s}"+".FINAL.vcf",s=samples),
            expand(config['project']['workpath']+"/vardict_out/{s}"+".FINAL.vcf",s=samples),
            expand(config['project']['workpath']+"/vardict_out/oncotator_out/{s}"+".maf",s=samples),
            expand(config['project']['workpath']+"/mutect_out/oncotator_out/{s}"+".maf",s=samples),
            expand(config['project']['workpath']+"/merged_somatic_variants/{s}"+".merged.vcf",s=samples),
            expand(config['project']['workpath']+"/merged_somatic_variants/oncotator_out/{s}"+".maf",s=samples),
            config['project']['workpath']+"/exome_targets.bed",
            expand("manta_out/{s}/results/variants/candidateSV.vcf.gz", s=samples),
            "sample_network_mqc.png",
            "admixture_out/admixture_table.tsv"
    output:
    params: rname="final"
    shell:  """
             module load multiqc/1.7; multiqc -f -e featureCounts .; mv *.out slurmfiles/; module load perl/5.18.4; perl Scripts/summarize_usage.pl; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped logfiles/

            """

elif config['project']['annotation'] == "hg38":

  rule all_exomeseq_somatic_tumoronly:
    input:  expand("{s}"+".recal.bam",s=samples),
            expand(config['project']['workpath']+"/mutect2_out/{s}"+".FINALmutect2.vcf",s=samples),
            expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{s}"+".maf",s=samples),
            config['project']['workpath']+"/mutect2_out/mutect2_maf_summary.pdf",
            config['project']['workpath']+"/mutect_out/mutect_maf_summary.pdf",
            config['project']['workpath']+"/vardict_out/vardict_maf_summary.pdf",
            config['project']['workpath']+"/merged_somatic_variants/merged_maf_summary.pdf",
            config['project']['workpath']+"/mutect2_out",
            "sample_network_mqc.png",
            config['project']['workpath']+"/mutect2_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/mutect_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/vardict_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/merged_somatic_variants/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/mutect_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/vardict_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/merged_somatic_variants/mutsigCV_out/somatic.sig_genes.txt",
            expand(config['project']['workpath']+"/mutect_out/{s}"+".FINAL.vcf",s=samples),
            expand(config['project']['workpath']+"/vardict_out/{s}"+".FINAL.vcf",s=samples),
            expand(config['project']['workpath']+"/vardict_out/oncotator_out/{s}"+".maf",s=samples),
            expand(config['project']['workpath']+"/mutect_out/oncotator_out/{s}"+".maf",s=samples),
            expand(config['project']['workpath']+"/merged_somatic_variants/{s}"+".merged.vcf",s=samples),
            config['project']['workpath']+"/exome_targets.bed",
            expand("manta_out/{s}/results/variants/candidateSV.vcf.gz", s=samples),
            "admixture_out/admixture_table.tsv"
    output:
    params: rname="final"
    shell:  """
             module load multiqc/1.7; multiqc -f -e featureCounts .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped logfiles/

            """
            
elif config['project']['annotation'] == "mm10":

  rule all_exomeseq_somatic_tumoronly:
    input:  expand("{s}"+".recal.bam",s=samples),
            expand(config['project']['workpath']+"/mutect2_out/{s}"+".FINALmutect2.vcf",s=samples),
            expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{s}"+".maf",s=samples),
            config['project']['workpath']+"/mutect2_out",
            config['project']['workpath']+"/mutect2_out/mutect2_maf_summary.pdf",
            config['project']['workpath']+"/mutect_out/mutect_maf_summary.pdf",
            config['project']['workpath']+"/vardict_out/vardict_maf_summary.pdf",
            config['project']['workpath']+"/merged_somatic_variants/merged_maf_summary.pdf",
            "sample_network_mqc.png",
            config['project']['workpath']+"/mutect2_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/mutect_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/vardict_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/merged_somatic_variants/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/mutect_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/vardict_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/merged_somatic_variants/mutsigCV_out/somatic.sig_genes.txt",
            expand(config['project']['workpath']+"/mutect_out/{s}"+".FINAL.vcf",s=samples),
            expand(config['project']['workpath']+"/vardict_out/{s}"+".FINAL.vcf",s=samples),
            expand(config['project']['workpath']+"/vardict_out/oncotator_out/{s}"+".maf",s=samples),
            expand(config['project']['workpath']+"/mutect_out/oncotator_out/{s}"+".maf",s=samples),
            expand(config['project']['workpath']+"/merged_somatic_variants/{s}"+".merged.vcf",s=samples),
            config['project']['workpath']+"/exome_targets.bed",
            expand("manta_out/{s}/results/variants/candidateSV.vcf.gz", s=samples),
            "admixture_out/admixture_table.tsv"
    output:
    params: rname="final"
    shell:  """
             module load multiqc/1.7; multiqc -f -e featureCounts .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped logfiles/

            """