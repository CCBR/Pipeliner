rule bwa_index_ref:
    input:  config['references']['GENOME']
    output: config['references']['GENOME']+".bwt"
    params: bwa=config['bin']['BWA'],ref=config['references']['GENOME'],rname="pl:bwaindex"
    shell:  "{params.bwa} index -a bwtsw {params.ref}"
