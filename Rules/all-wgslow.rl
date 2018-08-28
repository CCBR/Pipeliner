if config['project']['annotation'] == "hg19":
  rule all_wgslow:
    input: "combined.vcf",
           config['project']['workpath']+"/full_annot.txt.zip",
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
           config['project']['workpath']+"/svaba_out/svaba.log",
           expand("{s}"+".g.vcf",s=samples),
           expand("canvas_out/{s}/CNV.vcf.gz", s=samples),
           "admixture_out/admixture_table.tsv",

    output: 
    params: rname="final"
    shell:  """
             module load multiqc/1.4; multiqc -f -e featureCounts .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """

elif config['project']['annotation'] == "hg38":
  rule all_wgslow:
    input: "combined.vcf",
           config['project']['workpath']+"/full_annot.txt.zip",
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
           config['project']['workpath']+"/svaba_out/svaba.log",
           expand("{s}"+".g.vcf",s=samples),
           expand("canvas_out/{s}/CNV.vcf.gz", s=samples),
           "admixture_out/admixture_table.tsv"

    output: 
    params: rname="final"
    shell:  """
             module load multiqc/1.4; multiqc -f -e featureCounts .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """

elif config['project']['annotation'] == "mm10":
  rule all_wgslow:
    input: "combined.vcf",
           config['project']['workpath']+"/full_annot.txt.zip",
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
           config['project']['workpath']+"/svaba_out/svaba.log",
           expand("canvas_out/{s}/CNV.vcf.gz", s=samples),
           "admixture_out/admixture_table.tsv"
    output: 
    params: rname="final"
    shell:  """
             module load multiqc/1.4; multiqc -f -e featureCounts .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl; rm *realign.bai; mv distance.cluster0 distance.cluster1 distance.cluster2 distance.cluster3 distance.nosex samples.txt plink.map plink.ped *.avia_status.txt *.avia.log *_genotypes.vcf logfiles/

            """