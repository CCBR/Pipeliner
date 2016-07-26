rule index_ref:
    input:  config['references']['GENOME']
    output: config['references']['GENOME']+".bwt"
    params: bwa=config['bin']['BWA'],ref=config['references']['GENOME']
    shell:  "{params.bwa} index {params.ref}"
