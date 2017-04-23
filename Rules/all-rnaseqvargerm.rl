rule all_rnaseqvargerm:
    input: "combined.vcf",
           config['project']['workpath']+"/full_annot.txt.zip",
           "sample_network.bmp",
           "targets.strictFilter.snpeff.vcf",
           expand("sample_vcfs/{s}"+".stats.csv",s=samples),
           "targets.strictFilter.vcf",
           expand("QC/{s}.R1_fastqc.html",s=samples),
           expand("QC/{s}.{r}.trimmed_fastqc.html",s=samples,r=['R1','R1']),
           expand("QC/{s}.{r}.trimmed_screen.txt",s=samples,r=['R1','R1']),
           expand("QC/{s}.{r}.trimmed_screen.png",s=samples,r=['R1','R1']),
           expand("{s}.dedup.bam",s=samples),
#           expand("{s}.dedup.bam.onTarget.bam_stats",s=samples),
#           expand("{s}.dedup.bam.onTarget.bam",s=samples),
#           expand("{s}.dedup.bam.bam_stats",s=samples),
            expand("QC/{s}_run_trimmomatic.err",s=samples),
#            expand("QC/{s}.qualimapReport/genome_results.txt",s=samples),
#            config['project']['id']+"_"+config['project']['flowcellid']+".xlsx"
    output: 
    params: rname="final"
    shell:  """
             module load multiqc; multiqc -f -e featureCounts .; mv *.out slurmfiles/; mv *.fin.bam.intervals logfiles/; rm *realign.bai; rm *sorted.bam.bai; mv *sorted.txt logfiles/; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """