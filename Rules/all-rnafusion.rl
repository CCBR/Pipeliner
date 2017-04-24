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
    output: 
    params: rname="final"
    shell:  """
             module load multiqc; multiqc -f -e featureCounts .; mv *.out slurmfiles/

            """