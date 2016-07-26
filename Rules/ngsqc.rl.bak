rule ngsqc:
    input:  "{x}.R1."+config['project']['filetype'],
            "{x}.R2."+config['project']['filetype']
    output: "QC/{x}.R1."+config['project']['filetype']+"_filtered",
            "QC/{x}.R2."+config['project']['filetype']+"_filtered"
    params: illuqc=config['bin']['NGSQCTOOLKIT'],
            adapters=config['references']['adapter.file'],rname="pl:ngsqc"
    shell:  "{params.illuqc} -pe {input[0]} {input[1]} 3 A -o QC"
#    shell:  "{params.illuqc} -pe {input[0]} {input[1]} {params.adapters} A -o QC"



