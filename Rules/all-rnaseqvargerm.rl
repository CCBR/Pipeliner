rule all_rnaseqvargerm:
    input: "combined.vcf",
           "sample_network_mqc.png",
           "combined.strictFilter.snpeff.vcf",
           expand("sample_vcfs/{s}"+".stats.csv",s=samples),
           "combined.strictFilter.vcf",
           expand("QC/{s}.R1_fastqc.html",s=samples),
           expand("QC/{s}.{r}.trimmed_fastqc.html",s=samples,r=['R1','R1']),
           expand("QC/{s}.{r}.trimmed_screen.txt",s=samples,r=['R1','R1']),
           expand("QC/{s}.{r}.trimmed_screen.png",s=samples,r=['R1','R1']),
           expand("QC/{s}_run_trimmomatic.err",s=samples),
           "admixture_out/admixture_table.tsv"
    output: 
    params: rname="final"
    shell:  """
             module load multiqc/1.7; multiqc -f -e featureCounts .; mv *.out slurmfiles/; module load perl/5.18.4; perl Scripts/summarize_usage.pl; mv *.fin.bam.intervals logfiles/; rm *realign.bai; rm *sorted.bam.bai; mv *sorted.txt logfiles/; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """