rule all_exomeseq_germline:
    input: "combined.vcf.gz",
           "sample_network_mqc.png",
           "exome.snpeff.vcf",
           expand("sample_vcfs/{s}"+".stats.csv",s=samples),
           "exome_targets.bed",
           "exome.strictFilter.vcf.gz",
           "manta_out/results/variants/diploidSV.vcf.gz",
           "admixture_out/admixture_table.tsv"
    output: 
    params: rname="final"
    shell:  """
             module load multiqc/1.7; multiqc -f -e featureCounts .; mv *.out slurmfiles/; module load perl/5.18.4; perl Scripts/summarize_usage.pl; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """