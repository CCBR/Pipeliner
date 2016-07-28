rule cnvkit_somatic:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".realign.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".realign.bam"
    output: calls="cnvkit_out/{x}_calls.cns",
            gainloss="cnvkit_out/{x}_gainloss.tsv"
    params: tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1], normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0], targets=config['references'][pfamily]['CNVKIT_TARGETS'], antitargets=config['references'][pfamily]['CNVKIT_ANTITARGETS'], genome=config['references'][pfamily]['CNVKITGENOME'],rname="pl:cnvkit_somatic"
    shell: "module load cnvkit; cnvkit.py batch -n {input.normal} --fasta {params.genome} --targets {params.targets} --antitargets {params.antitargets} --output-reference cnvkit_out/{params.normalsample}.reference.cnn --output-dir cnvkit_out {input.tumor}; cnvkit.py call -o {output.calls} cnvkit_out/{params.tumorsample}.cnr; cnvkit.py gainloss -s cnvkit_out/{params.tumorsample}.cns --drop-low-coverage -o {output.gainloss} cnvkit_out/{params.tumorsample}.cnr"