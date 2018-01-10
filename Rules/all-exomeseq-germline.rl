rule all_exomeseq_germline:
    input: "combined.vcf",
           "exome.recode.vcf",
           config['project']['workpath']+"/full_annot.txt.zip",
#           "variants.database",
           "sample_network.bmp",
           "exome.snpeff.vcf",
           expand("sample_vcfs/{s}"+".stats.csv",s=samples),
           "exome_targets.bed",
           "exome.strictFilter.vcf",
           "manta_out/results/variants/diploidSV.vcf.gz",
           "admixture_out/admixture_table.tsv"
    output: 
    params: rname="final"
    shell:  """
             module load multiqc/1.3; multiqc -f -e featureCounts .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """