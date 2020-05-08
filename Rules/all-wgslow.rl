if config['project']['annotation'] == "hg19":
  rule all_wgslow:
    input: "combined.vcf",
           "sample_network_mqc.png",
           expand("sample_vcfs/{s}"+".stats.csv",s=samples),           
           "combined.snpeff.vcf",
           "combined.strictFilter.vcf",
           "manta_out/results/variants/diploidSV.vcf.gz",
           config['project']['workpath']+"/svaba_out/svaba.log",
#           expand("{s}"+".g.vcf",s=samples),
           expand("canvas_out/{s}/CNV.vcf.gz", s=samples),
           "admixture_out/admixture_table.tsv",
           expand("{s}.recal.bam",s=samples),
           expand("QC/{s}.R1_fastqc.html",s=samples),
#          expand("QC/{s}.R1."+config['project']['filetype']+"_filtered",s=samples),
#          expand("{s}.R1.trimmed.fastq.bz2",s=samples),
#          expand("QC/{s}.qualimapReport",s=samples),
           expand("QC/{s}.{r}.trimmed_fastqc.html",s=samples,r=['R1','R1']),
           expand("QC/{s}.{r}.trimmed_screen.txt",s=samples,r=['R1','R1']),
           expand("QC/{s}.{r}.trimmed_screen.png",s=samples,r=['R1','R1']),
           expand("QC/{s}_run_trimmomatic.err",s=samples),
           expand("QC/{s}.qualimapReport/genome_results.txt",s=samples),

    output: 
    params: rname="final"
    shell:  """
             module load multiqc/1.7; multiqc -f -e featureCounts .; mv *.out slurmfiles/; module load perl/5.18.4; perl Scripts/summarize_usage.pl; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """

elif config['project']['annotation'] == "hg38":
  rule all_wgslow:
    input: "combined.vcf",
           "sample_network_mqc.png",
           expand("sample_vcfs/{s}"+".stats.csv",s=samples),           
           "combined.snpeff.vcf",
           "combined.strictFilter.vcf",
           "manta_out/results/variants/diploidSV.vcf.gz",
           config['project']['workpath']+"/svaba_out/svaba.log",
#           expand("{s}"+".g.vcf",s=samples),
           expand("canvas_out/{s}/CNV.vcf.gz", s=samples),
           "admixture_out/admixture_table.tsv",
           expand("{s}.recal.bam",s=samples),
           expand("QC/{s}.R1_fastqc.html",s=samples),
#          expand("QC/{s}.R1."+config['project']['filetype']+"_filtered",s=samples),
#          expand("{s}.R1.trimmed.fastq.bz2",s=samples),
#          expand("QC/{s}.qualimapReport",s=samples),
           expand("QC/{s}.{r}.trimmed_fastqc.html",s=samples,r=['R1','R1']),
           expand("QC/{s}.{r}.trimmed_screen.txt",s=samples,r=['R1','R1']),
           expand("QC/{s}.{r}.trimmed_screen.png",s=samples,r=['R1','R1']),
           expand("QC/{s}_run_trimmomatic.err",s=samples),
           expand("QC/{s}.qualimapReport/genome_results.txt",s=samples),

    output: 
    params: rname="final"
    shell:  """
             module load multiqc/1.7; multiqc -f -e featureCounts .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """

elif config['project']['annotation'] == "mm10":
  rule all_wgslow:
    input: "combined.vcf",
           "sample_network_mqc.png",
           expand("sample_vcfs/{s}"+".stats.csv",s=samples),           
           "combined.snpeff.vcf",
           "combined.strictFilter.vcf",
           "manta_out/results/variants/diploidSV.vcf.gz",
           config['project']['workpath']+"/svaba_out/svaba.log",
           expand("canvas_out/{s}/CNV.vcf.gz", s=samples),
           "admixture_out/admixture_table.tsv",
           expand("{s}.recal.bam",s=samples),
           expand("QC/{s}.R1_fastqc.html",s=samples),
#          expand("QC/{s}.R1."+config['project']['filetype']+"_filtered",s=samples),
#          expand("{s}.R1.trimmed.fastq.bz2",s=samples),
#          expand("QC/{s}.qualimapReport",s=samples),
           expand("QC/{s}.{r}.trimmed_fastqc.html",s=samples,r=['R1','R1']),
           expand("QC/{s}.{r}.trimmed_screen.txt",s=samples,r=['R1','R1']),
           expand("QC/{s}.{r}.trimmed_screen.png",s=samples,r=['R1','R1']),
           expand("QC/{s}_run_trimmomatic.err",s=samples),
           expand("QC/{s}.qualimapReport/genome_results.txt",s=samples),
    output: 
    params: rname="final"
    shell:  """
             module load multiqc/1.7; multiqc -f -e featureCounts .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """