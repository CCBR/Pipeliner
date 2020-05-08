if config['project']['annotation'] == "hg19":
  rule all_rnafusion:
    input: expand("starfusion/oncofuse/{x}/{x}.oncofuse.input",x=samples),
           expand("starfusion/oncofuse/{x}/{x}.oncofuse.output",x=samples),
           expand("fusioncatcher/oncofuse/{x}/{x}.oncofuse.input",x=samples),
           expand("fusioncatcher/oncofuse/{x}/{x}.oncofuse.output",x=samples),
           expand("QC/{x}.R1.trimmed_fastqc.html",x=samples),
           expand("QC/{x}.R2.trimmed_fastqc.html",x=samples),
           expand("QC/{x}.R1.trimmed_screen.txt",x=samples),
           expand("QC/{x}.R1.trimmed_screen.png",x=samples),
           expand("QC/{x}.R2.trimmed_screen.txt",x=samples),
           expand("QC/{x}.R2.trimmed_screen.png",x=samples),
           expand("QC/{x}.R1_fastqc.html",x=samples),
           expand("QC/{x}.R2_fastqc.html",x=samples),
           config['project']['workpath']+"/fusioncatcher",
           config['project']['workpath']+"/fusioncatcher/fusioninspector",
           config['project']['workpath']+"/fusioncatcher/oncofuse",
           config['project']['workpath']+"/starfusion",
           config['project']['workpath']+"/starfusion/fusioninspector",
           config['project']['workpath']+"/starfusion/oncofuse",
           expand("fusioncatcher/{x}/final-list_candidate-fusion-genes.txt",x=samples),
           expand("starfusion/fusioninspector/{x}/{x}.fusion_predictions.final",x=samples),
           expand("fusioncatcher/fusioninspector/{x}/{x}.fusion_predictions.final",x=samples),
           expand("starfusion/{x}/star-fusion.fusion_predictions.tsv",x=samples),
#           expand("{name}.RnaSeqMetrics.txt",name=samples),
#           "expression/RawCountFile_genes_filtered.txt",
#           "expression/RawCountFile_junctions_filtered.txt",
#           "expression/RawCountFile_genejunctions_filtered.txt",
#           expand("expression/{name}.star.count.overlap.txt",name=samples),
#           "expression/RawCountFileOverlap.txt",
#           "expression/RawCountFileStar.txt",
#           expand("{name}.rsem.genes.results",name=samples),
#           "RawCountFile_RSEM_genes_filtered.txt",
#           expand("QC/{x}_readlength.txt",x=samples),
           expand("fusioninspector/{x}/{x}.fusion_predictions.final",x=samples),
           expand("oncofuse/{x}/{x}.oncofuse.output",x=samples)
    output: 
    params: rname="final"
    shell:  """
             Scripts/fusionSummary.sh; module load multiqc/1.7; multiqc -f .; rm *featureCounts; mv *.out slurmfiles/; module load perl/5.18.4; perl Scripts/summarize_usage.pl

            """

elif config['project']['annotation'] == "hg38":
  rule all_rnafusion:
    input: expand("starfusion/oncofuse/{x}/{x}.oncofuse.input",x=samples),
           expand("starfusion/oncofuse/{x}/{x}.oncofuse.output",x=samples),
           expand("fusioncatcher/oncofuse/{x}/{x}.oncofuse.input",x=samples),
           expand("fusioncatcher/oncofuse/{x}/{x}.oncofuse.output",x=samples),
           expand("QC/{x}.R1.trimmed_fastqc.html",x=samples),
           expand("QC/{x}.R2.trimmed_fastqc.html",x=samples),
           expand("QC/{x}.R1.trimmed_screen.txt",x=samples),
           expand("QC/{x}.R1.trimmed_screen.png",x=samples),
           expand("QC/{x}.R2.trimmed_screen.txt",x=samples),
           expand("QC/{x}.R2.trimmed_screen.png",x=samples),
           expand("QC/{x}.R1_fastqc.html",x=samples),
           expand("QC/{x}.R2_fastqc.html",x=samples),
           config['project']['workpath']+"/fusioncatcher",
           config['project']['workpath']+"/fusioncatcher/fusioninspector",
           config['project']['workpath']+"/fusioncatcher/oncofuse",
           config['project']['workpath']+"/starfusion",
           config['project']['workpath']+"/starfusion/fusioninspector",
           config['project']['workpath']+"/starfusion/oncofuse",
           expand("fusioncatcher/{x}/final-list_candidate-fusion-genes.txt",x=samples),
           expand("starfusion/fusioninspector/{x}/{x}.fusion_predictions.final",x=samples),
           expand("fusioncatcher/fusioninspector/{x}/{x}.fusion_predictions.final",x=samples),
           expand("starfusion/{x}/star-fusion.fusion_predictions.tsv",x=samples),
#           expand("{name}.RnaSeqMetrics.txt",name=samples),
#           "expression/RawCountFile_genes_filtered.txt",
#           "expression/RawCountFile_junctions_filtered.txt",
#           "expression/RawCountFile_genejunctions_filtered.txt",
#           expand("expression/{name}.star.count.overlap.txt",name=samples),
#           "expression/RawCountFileOverlap.txt",
#           "expression/RawCountFileStar.txt",
#           expand("{name}.rsem.genes.results",name=samples),
#           "RawCountFile_RSEM_genes_filtered.txt",
#           expand("QC/{x}_readlength.txt",x=samples),
           expand("fusioninspector/{x}/{x}.fusion_predictions.final",x=samples),
           expand("oncofuse/{x}/{x}.oncofuse.output",x=samples)
    output: 
    params: rname="final"
    shell:  """
             Scripts/fusionSummary.sh; module load multiqc/1.7; multiqc -f .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl

            """

elif config['project']['annotation'] == "mm10":
  rule all_rnafusion:
    input: expand("QC/{x}.R1.trimmed_fastqc.html",x=samples),
           expand("QC/{x}.R2.trimmed_fastqc.html",x=samples),
           expand("QC/{x}.R1.trimmed_screen.txt",x=samples),
           expand("QC/{x}.R1.trimmed_screen.png",x=samples),
           expand("QC/{x}.R2.trimmed_screen.txt",x=samples),
           expand("QC/{x}.R2.trimmed_screen.png",x=samples),
           expand("QC/{x}.R1_fastqc.html",x=samples),
           expand("QC/{x}.R2_fastqc.html",x=samples),
           config['project']['workpath']+"/fusioncatcher",
           config['project']['workpath']+"/fusioncatcher/fusioninspector",
           config['project']['workpath']+"/fusioncatcher/oncofuse",
           config['project']['workpath']+"/starfusion",
           config['project']['workpath']+"/starfusion/fusioninspector",
           config['project']['workpath']+"/starfusion/oncofuse",
           expand("fusioncatcher/{x}/final-list_candidate-fusion-genes.txt",x=samples),
           expand("starfusion/fusioninspector/{x}/{x}.fusion_predictions.final",x=samples),
           expand("fusioncatcher/fusioninspector/{x}/{x}.fusion_predictions.final",x=samples),
           expand("starfusion/{x}/star-fusion.fusion_predictions.tsv",x=samples),
#           expand("{name}.RnaSeqMetrics.txt",name=samples),
#           "expression/RawCountFile_genes_filtered.txt",
#           "expression/RawCountFile_junctions_filtered.txt",
#           "expression/RawCountFile_genejunctions_filtered.txt",
#           expand("expression/{name}.star.count.overlap.txt",name=samples),
#           "expression/RawCountFileOverlap.txt",
#           "expression/RawCountFileStar.txt",
#           expand("{name}.rsem.genes.results",name=samples),
#           "RawCountFile_RSEM_genes_filtered.txt",
#           expand("QC/{x}_readlength.txt",x=samples),
           expand("fusioninspector/{x}/{x}.fusion_predictions.final",x=samples),
           expand("oncofuse/{x}/{x}.oncofuse.output",x=samples)
    output: 
    params: rname="final"
    shell:  """
             Scripts/fusionSummary.sh; module load multiqc/1.7; multiqc -f .; mv *.out slurmfiles/; perl Scripts/summarize_usage.pl

            """