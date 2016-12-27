rule all_exomeseq_somatic:
    input: "combined.vcf",
           "full_annot.txt.zip",
           "variants.database",
           "sample_network.bmp",
           "targets.snpeff.vcf",
           expand("sample_vcfs/{s}"+".stats.csv",s=samples),
           "targets.strictFilter.vcf"
    output: 
    params: rname="final"
    shell:  """
             module load multiqc; multiqc -f -e featureCounts .; mv *.out slurmfiles/; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """