if config['project']['annotation'] == "hg19":

  rule all_exomeseq_somatic:
    input:  expand("{s}"+".recal.bam",s=samples),
            expand(config['project']['workpath']+"/mutect_out/{p}"+".FINAL.vcf",p=pairs),
            expand(config['project']['workpath']+"/mutect2_out/{p}"+".FINALmutect2.vcf",p=pairs),
            expand(config['project']['workpath']+"/strelka_out/{p}"+"_FINAL.vcf",p=pairs),
            expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/strelka_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/mutect_out/oncotator_out/{p}"+".maf",p=pairs),
            config['project']['workpath']+"/mutect2_out/mutect2_maf_summary.pdf",
            config['project']['workpath']+"/mutect_out/mutect_maf_summary.pdf",
            config['project']['workpath']+"/strelka_out/strelka_maf_summary.pdf",
            expand(config['project']['workpath']+"/conpair_out/{p}.conpair", p=pairs),
            config['project']['workpath']+"/mutect2_out/mutect2_variants.database",
            config['project']['workpath']+"/mutect_out/mutect_variants.database",
            config['project']['workpath']+"/strelka_out/strelka_variants.database",
            config['project']['workpath']+"/mutect2_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/strelka_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/strelka_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/mutect_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/mutect_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/exome_targets.bed",
            expand("manta_out/{p}/results/variants/candidateSV.vcf.gz", p=pairs),
            "admixture_out/admixture_table.tsv",
            "QC/decoy",
            expand("sequenza_out/{p}"+"_alternative_solutions.txt",p=pairs),
            expand("freec_out/pass2/{p}"+".recal.bam_CNVs.p.value.txt",p=pairs),
            "sample_network_mqc.png",
            expand(config['project']['workpath']+"/vardict_out/{p}"+".FINAL.vcf",p=pairs),
            expand(config['project']['workpath']+"/vardict_out/oncotator_out/{p}"+".maf",p=pairs),
            config['project']['workpath']+"/vardict_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/merged_somatic_variants/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/vardict_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/merged_somatic_variants/mutsigCV_out/somatic.sig_genes.txt",
            expand(config['project']['workpath']+"/merged_somatic_variants/{p}"+".merged.vcf",p=pairs),
            expand(config['project']['workpath']+"/merged_somatic_variants/oncotator_out/{p}"+".maf",p=pairs),
            config['project']['workpath']+"/vardict_out/vardict_maf_summary.pdf",
            config['project']['workpath']+"/merged_somatic_variants/merged_maf_summary.pdf",
    output:
    params: rname="final"
    shell:  """
             module load multiqc/1.7; multiqc -f -e featureCounts .; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped logfiles/

            """
elif config['project']['annotation'] == "hg38":

  rule all_exomeseq_somatic:
    input:  expand("{s}"+".recal.bam",s=samples),
            expand(config['project']['workpath']+"/mutect_out/{p}"+".FINAL.vcf",p=pairs),
            expand(config['project']['workpath']+"/mutect2_out/{p}"+".FINALmutect2.vcf",p=pairs),
            expand(config['project']['workpath']+"/strelka_out/{p}"+"_FINAL.vcf",p=pairs),
            expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/strelka_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/mutect_out/oncotator_out/{p}"+".maf",p=pairs),
            config['project']['workpath']+"/mutect2_out/mutect2_maf_summary.pdf",
            config['project']['workpath']+"/mutect_out/mutect_maf_summary.pdf",
            config['project']['workpath']+"/strelka_out/strelka_maf_summary.pdf",
            expand(config['project']['workpath']+"/conpair_out/{p}.conpair", p=pairs),
            config['project']['workpath']+"/mutect2_out/mutect2_variants.database",
            config['project']['workpath']+"/mutect_out/mutect_variants.database",
            config['project']['workpath']+"/strelka_out/strelka_variants.database",
            config['project']['workpath']+"/mutect2_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/strelka_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/strelka_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/mutect_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/mutect_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/exome_targets.bed",
            expand("manta_out/{p}/results/variants/candidateSV.vcf.gz", p=pairs),
            "admixture_out/admixture_table.tsv",
            "QC/decoy",
            expand("sequenza_out/{p}"+"_alternative_solutions.txt",p=pairs),
            expand("freec_out/pass2/{p}"+".recal.bam_CNVs.p.value.txt",p=pairs),
            "sample_network_mqc.png",
            expand(config['project']['workpath']+"/vardict_out/{p}"+".FINAL.vcf",p=pairs),
            expand(config['project']['workpath']+"/vardict_out/oncotator_out/{p}"+".maf",p=pairs),
            config['project']['workpath']+"/vardict_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/merged_somatic_variants/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/vardict_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/merged_somatic_variants/mutsigCV_out/somatic.sig_genes.txt",
            expand(config['project']['workpath']+"/merged_somatic_variants/{p}"+".merged.vcf",p=pairs),
            expand(config['project']['workpath']+"/merged_somatic_variants/oncotator_out/{p}"+".maf",p=pairs),
            config['project']['workpath']+"/vardict_out/vardict_maf_summary.pdf",
            config['project']['workpath']+"/merged_somatic_variants/merged_maf_summary.pdf",
    output:
    params: rname="final"
    shell:  """
             module load multiqc/1.7; multiqc -f -e featureCounts .; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped logfiles/

            """

elif config['project']['annotation'] == "mm10":

  rule all_exomeseq_somatic:
    input:  expand("{s}"+".recal.bam",s=samples),
            expand(config['project']['workpath']+"/mutect_out/{p}"+".FINAL.vcf",p=pairs),
            expand(config['project']['workpath']+"/mutect2_out/{p}"+".FINALmutect2.vcf",p=pairs),
            expand(config['project']['workpath']+"/strelka_out/{p}"+"_FINAL.vcf",p=pairs),
            expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/strelka_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/mutect_out/oncotator_out/{p}"+".maf",p=pairs),
            expand(config['project']['workpath']+"/conpair_out/{p}.conpair", p=pairs),
            config['project']['workpath']+"/mutect2_out/mutect2_variants.database",
            config['project']['workpath']+"/mutect_out/mutect_variants.database",
            config['project']['workpath']+"/strelka_out/strelka_variants.database",
            config['project']['workpath']+"/mutect2_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/strelka_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/strelka_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/mutect_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/mutect_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/exome_targets.bed",
            expand("manta_out/{p}/results/variants/candidateSV.vcf.gz", p=pairs),
            "admixture_out/admixture_table.tsv",
            "QC/decoy",
            "sample_network_mqc.png",
            expand(config['project']['workpath']+"/vardict_out/{p}"+".FINAL.vcf",p=pairs),
            expand(config['project']['workpath']+"/vardict_out/oncotator_out/{p}"+".maf",p=pairs),
            config['project']['workpath']+"/vardict_out/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/merged_somatic_variants/oncotator_out/final_filtered.maf",
            config['project']['workpath']+"/vardict_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/merged_somatic_variants/mutsigCV_out/somatic.sig_genes.txt",
            expand(config['project']['workpath']+"/merged_somatic_variants/{p}"+".merged.vcf",p=pairs),
            expand(config['project']['workpath']+"/merged_somatic_variants/oncotator_out/{p}"+".maf",p=pairs),
            config['project']['workpath']+"/vardict_out/vardict_maf_summary.pdf",
            config['project']['workpath']+"/merged_somatic_variants/merged_maf_summary.pdf",
    output:
    params: rname="final"
    shell:  """
             module load multiqc/1.7; multiqc -f -e featureCounts .; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped logfiles/

            """
