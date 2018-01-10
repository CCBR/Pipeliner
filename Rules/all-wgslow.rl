if config['project']['annotation'] == "hg19":
  rule all_wgslow:
    input: "combined.vcf",
           config['project']['workpath']+"/full_annot.txt.zip",
#           "variants.database",
           "sample_network.bmp",
           expand("sample_vcfs/{s}"+".stats.csv",s=samples),           
           "combined.snpeff.vcf",
           "combined.strictFilter.vcf",
           "manta_out/results/variants/diploidSV.vcf.gz",
           config['project']['workpath']+"/delly_out/deletions.bcf",
           config['project']['workpath']+"/delly_out/inversions.bcf",
           config['project']['workpath']+"/delly_out/translocations.bcf",
           config['project']['workpath']+"/delly_out/duplications.bcf",
           config['project']['workpath']+"/delly_out/insertions.bcf",
#           "pindel_out/pindel_calls_chr1_INT_final",
#           "pindel_out/pindel_calls_chr2_INT_final",
#           "pindel_out/pindel_calls_chr3_INT_final",
#           "pindel_out/pindel_calls_chr4_INT_final",
#           "pindel_out/pindel_calls_chr5_INT_final",
#           "pindel_out/pindel_calls_chr6_INT_final",
#           "pindel_out/pindel_calls_chr7_INT_final",
#           "pindel_out/pindel_calls_chr8_INT_final",
#           "pindel_out/pindel_calls_chr9_INT_final",
#           "pindel_out/pindel_calls_chr10_INT_final",
#           "pindel_out/pindel_calls_chr11_INT_final",
#           "pindel_out/pindel_calls_chr12_INT_final",
#           "pindel_out/pindel_calls_chr13_INT_final",
#           "pindel_out/pindel_calls_chr14_INT_final",
#           "pindel_out/pindel_calls_chr15_INT_final",
#           "pindel_out/pindel_calls_chr16_INT_final",
#           "pindel_out/pindel_calls_chr17_INT_final",
#           "pindel_out/pindel_calls_chr18_INT_final",
#           "pindel_out/pindel_calls_chr19_INT_final",
#           "pindel_out/pindel_calls_chr20_INT_final",
#           "pindel_out/pindel_calls_chr21_INT_final",
#           "pindel_out/pindel_calls_chr22_INT_final",
#           "pindel_out/pindel_calls_chrX_INT_final",
#           "pindel_out/pindel_calls_chrY_INT_final",
#           "cnvkit_out/cnvkit_heatmap.pdf",
#           "breakdancer_out/file.ctx",
           config['project']['workpath']+"/svaba_out/svaba.log",
           expand("{s}"+".g.vcf",s=samples),
           expand("canvas_out/{s}/CNV.vcf.gz", s=samples),
           "admixture_out/admixture_table.tsv",

    output: 
    params: rname="final"
    shell:  """
             module load multiqc/1.3; multiqc -f -e featureCounts .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """

elif config['project']['annotation'] == "hg38":
  rule all_wgslow:
    input: "combined.vcf",
           config['project']['workpath']+"/full_annot.txt.zip",
#           "variants.database",
           "sample_network.bmp",
           expand("sample_vcfs/{s}"+".stats.csv",s=samples),           
           "combined.snpeff.vcf",
           "combined.strictFilter.vcf",
           "manta_out/results/variants/diploidSV.vcf.gz",
           config['project']['workpath']+"/delly_out/deletions.bcf",
           config['project']['workpath']+"/delly_out/inversions.bcf",
           config['project']['workpath']+"/delly_out/translocations.bcf",
           config['project']['workpath']+"/delly_out/duplications.bcf",
           config['project']['workpath']+"/delly_out/insertions.bcf",
#           "pindel_out/pindel_calls_chr1_INT_final",
#           "pindel_out/pindel_calls_chr2_INT_final",
#           "pindel_out/pindel_calls_chr3_INT_final",
#           "pindel_out/pindel_calls_chr4_INT_final",
#           "pindel_out/pindel_calls_chr5_INT_final",
#           "pindel_out/pindel_calls_chr6_INT_final",
#           "pindel_out/pindel_calls_chr7_INT_final",
#           "pindel_out/pindel_calls_chr8_INT_final",
#           "pindel_out/pindel_calls_chr9_INT_final",
#           "pindel_out/pindel_calls_chr10_INT_final",
#           "pindel_out/pindel_calls_chr11_INT_final",
#           "pindel_out/pindel_calls_chr12_INT_final",
#           "pindel_out/pindel_calls_chr13_INT_final",
#           "pindel_out/pindel_calls_chr14_INT_final",
#           "pindel_out/pindel_calls_chr15_INT_final",
#           "pindel_out/pindel_calls_chr16_INT_final",
#           "pindel_out/pindel_calls_chr17_INT_final",
#           "pindel_out/pindel_calls_chr18_INT_final",
#           "pindel_out/pindel_calls_chr19_INT_final",
#           "pindel_out/pindel_calls_chr20_INT_final",
#           "pindel_out/pindel_calls_chr21_INT_final",
#           "pindel_out/pindel_calls_chr22_INT_final",
#           "pindel_out/pindel_calls_chrX_INT_final",
#           "pindel_out/pindel_calls_chrY_INT_final",
#           "cnvkit_out/cnvkit_heatmap.pdf",
#           "breakdancer_out/file.ctx",
           config['project']['workpath']+"/svaba_out/svaba.log",
           expand("{s}"+".g.vcf",s=samples),
           expand("canvas_out/{s}/CNV.vcf.gz", s=samples),
           "admixture_out/admixture_table.tsv"

    output: 
    params: rname="final"
    shell:  """
             module load multiqc/1.3; multiqc -f -e featureCounts .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """

elif config['project']['annotation'] == "mm10":
  rule all_wgslow:
    input: "combined.vcf",
           config['project']['workpath']+"/full_annot.txt.zip",
#           "variants.database",
           "sample_network.bmp",
           expand("sample_vcfs/{s}"+".stats.csv",s=samples),           
           "combined.snpeff.vcf",
           "combined.strictFilter.vcf",
           "manta_out/results/variants/diploidSV.vcf.gz",
           config['project']['workpath']+"/delly_out/deletions.bcf",
           config['project']['workpath']+"/delly_out/inversions.bcf",
           config['project']['workpath']+"/delly_out/translocations.bcf",
           config['project']['workpath']+"/delly_out/duplications.bcf",
           config['project']['workpath']+"/delly_out/insertions.bcf",
#           "pindel_out/pindel_calls_chr1_INT_final",
#           "pindel_out/pindel_calls_chr2_INT_final",
#           "pindel_out/pindel_calls_chr3_INT_final",
#           "pindel_out/pindel_calls_chr4_INT_final",
#           "pindel_out/pindel_calls_chr5_INT_final",
#           "pindel_out/pindel_calls_chr6_INT_final",
#           "pindel_out/pindel_calls_chr7_INT_final",
#           "pindel_out/pindel_calls_chr8_INT_final",
#           "pindel_out/pindel_calls_chr9_INT_final",
#           "pindel_out/pindel_calls_chr10_INT_final",
#           "pindel_out/pindel_calls_chr11_INT_final",
#           "pindel_out/pindel_calls_chr12_INT_final",
#           "pindel_out/pindel_calls_chr13_INT_final",
#           "pindel_out/pindel_calls_chr14_INT_final",
#           "pindel_out/pindel_calls_chr15_INT_final",
#           "pindel_out/pindel_calls_chr16_INT_final",
#           "pindel_out/pindel_calls_chr17_INT_final",
#           "pindel_out/pindel_calls_chr18_INT_final",
#           "pindel_out/pindel_calls_chr19_INT_final",
#           "pindel_out/pindel_calls_chrX_INT_final",
#           "pindel_out/pindel_calls_chrY_INT_final",
#           "cnvkit_out/cnvkit_heatmap.pdf",
#           "breakdancer_out/file.ctx",
           config['project']['workpath']+"/svaba_out/svaba.log",
           expand("canvas_out/{s}/CNV.vcf.gz", s=samples),
           "admixture_out/admixture_table.tsv"
    output: 
    params: rname="final"
    shell:  """
             module load multiqc/1.3; multiqc -f -e featureCounts .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """