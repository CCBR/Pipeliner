rule cnvkit_germ:
    input: lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam"
    output: cnr="cnvkit/"+lambda wildcards: config['project']['units'][wildcards.x]+".cnr",
            cns="cnvkit/"+lambda wildcards: config['project']['units'][wildcards.x]+".cns",
            heatmap="cnvkit/germline_cnvkit.heatmap"
    params: targets="cnvkit_targets.bed",antitargets="cnvkit_antitargets.bed",genome=config['references'][pfamily]['GENOME'],rname="pl:cnvkit_germ"
    shell: "module load cnvkit; cnvkit.py batch *.realign.bam -n -c --split --fasta {params.genome} --targets {params.targets} --antitargets {params.antitargets} --output-reference cnvkit/flatref.cnn --output-dir cnvkit/ --diagram; cnvkit.py heatmap cnvkit/*.cns -d -o cnvkit/germline_cnvkit.heatmap"