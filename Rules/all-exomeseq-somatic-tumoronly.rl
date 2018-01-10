if config['project']['annotation'] == "hg19":

  rule all_exomeseq_somatic_tumoronly:
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
            config['project']['workpath']+"/cnvkit_out/CNVkit_summary_heatmap.pdf",
            config['project']['workpath']+"/sample_network.bmp",
            config['project']['workpath']+"/mutect2_out/mutect2_variants.database",
            config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_filtered.maf",
            config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/exome_targets.bed",
            expand("manta_out/{s}/results/variants/candidateSV.vcf.gz", s=samples),
            "admixture_out/admixture_table.tsv"
    output:
    params: rname="final"
    shell:  """
             module load multiqc/1.3; multiqc -f -e featureCounts .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped logfiles/

            """

elif config['project']['annotation'] == "hg38":

  rule all_exomeseq_somatic_tumoronly:
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
            config['project']['workpath']+"/cnvkit_out/CNVkit_summary_heatmap.pdf",
            config['project']['workpath']+"/sample_network.bmp",
            config['project']['workpath']+"/mutect2_out/mutect2_variants.database",
            config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_filtered.maf",
            config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/exome_targets.bed",
            expand("manta_out/{s}/results/variants/candidateSV.vcf.gz", s=samples),
            "admixture_out/admixture_table.tsv"
    output:
    params: rname="final"
    shell:  """
             module load multiqc/1.3; multiqc -f -e featureCounts .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped logfiles/

            """
            
elif config['project']['annotation'] == "mm10":

  rule all_exomeseq_somatic_tumoronly:
    input:  expand("{s}"+".recal.bam",s=samples),
            expand(config['project']['workpath']+"/mutect2_out/{s}"+".FINALmutect2.vcf",s=samples),
            expand(config['project']['workpath']+"/mutect2_out/oncotator_out/{s}"+".maf",s=samples),
            config['project']['workpath']+"/mutect2_out",
            config['project']['workpath']+"/cnvkit_out",
            config['project']['workpath']+"/delly_out",
            expand(config['project']['workpath']+"/cnvkit_out/{s}_calls.cns", s=samples),
            expand(config['project']['workpath']+"/cnvkit_out/{s}_gainloss.tsv", s=samples),                        
            config['project']['workpath']+"/cnvkit_out/CNVkit_summary_heatmap.pdf",
            config['project']['workpath']+"/mutect2_out/mutect2_oncoplot.pdf",
            config['project']['workpath']+"/mutect2_out/mutect2_variants.database",
            config['project']['workpath']+"/sample_network.bmp",
#            config['project']['workpath']+"/mutect2_out/oncotator_out/mutect2_filtered.maf",
            config['project']['workpath']+"/mutect2_out/mutsigCV_out/somatic.sig_genes.txt",
            config['project']['workpath']+"/exome_targets.bed",
            expand("manta_out/{s}/results/variants/candidateSV.vcf.gz", s=samples),
            "admixture_out/admixture_table.tsv"
    output:
    params: rname="final"
    shell:  """
             module load multiqc/1.3; multiqc -f -e featureCounts .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped logfiles/

            """