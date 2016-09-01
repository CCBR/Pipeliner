rule all_exomeseq_germline:
    input: "combined.vcf",
           "exome.recode.vcf",
           "full_annot.txt.zip",
           "variants.database",
           "sample_network.bmp",
           "exome.snpeff.vcf"
    output: 
    params: rname="final"
    shell:  """
             mv *.out slurmfiles/; mv *.fin.bam.intervals logfiles/; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped viz-ccbr-cvdz5405.avia_status.txt viz-ccbr-cvdz5405.avia.log exome_genotypes.vcf logfiles/

            """