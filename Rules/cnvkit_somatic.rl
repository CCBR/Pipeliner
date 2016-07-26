rule cnvkit_somatic:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".realign.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".realign.bam"
    output: cnr="cnvkit_out/{x}.cnr",
            cns="cnvkit_out/{x}.cns"
    params: tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1], normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0], targets=config['references'][pfamily]['CNVKIT_TARGETS'], antitargets=config['references'][pfamily]['CNVKIT_ANTITARGETS'], genome=config['references'][pfamily]['CNVKITGENOME'],rname="pl:cnvkit_somatic"
    shell: "module load cnvkit; cnvkit.py batch -n {input.normal} --fasta {params.genome} --targets {params.targets} --antitargets {params.antitargets} --output-reference cnvkit_out/{params.normalsample}.reference.cnn --output-dir cnvkit_out {input.tumor}; mv cnvkit_out/{params.tumorsample}.cnr {output.cnr}; mv cnvkit_out/{params.normalsample}.cns {output.cns}"