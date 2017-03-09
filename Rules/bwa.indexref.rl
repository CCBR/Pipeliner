rule bwa_index_ref:
    input:  config['references'][pfamily]['GENOME']
    output: config['references'][pfamily]['GENOME']+".bwt"
    params: bwa=config['bin'][pfamily]['BWA'],ref=config['references'][pfamily]['GENOME'],rname="pl:bwaindex"
    shell:  "{params.bwa} index -a bwtsw {params.ref}"
