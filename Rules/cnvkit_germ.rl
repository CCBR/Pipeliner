rule cnvkit_germ:
    input: bams=expand(lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam")
    output: cnr=expand("cnvkit/"+lambda wildcards: config['project']['units'][wildcards.x]+".cnr"),
            cns=expand("cnvkit/"+lambda wildcards: config['project']['units'][wildcards.x]+".cns"),
            heatmap="cnvkit/germline_cnvkit.heatmap"
    params: targets="cnvkit_targets.bed",antitargets="cnvkit_antitargets.bed",genome=config['references'][pfamily]['GENOME'],rname="pl:cnvkit_germ"
    shell: "module load cnvkit/0.8.5; cnvkit.py batch --method wgs -r hs37d5_flatRef.cnn *.recal.bam; cnvkit.py heatmap -d -o cnvkit_heatmap.pdf *.cns