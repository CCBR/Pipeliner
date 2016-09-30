rule cnvkit_somatic:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="cnvkit_targets.bed"
    output: calls="cnvkit_out/{x}_calls.cns",
            gainloss="cnvkit_out/{x}_gainloss.tsv",
            dir="cnvkit_out/{x}_cnvkit"
    params: tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],access=config['references'][pfamily]['CNVKITACCESS'],targets="cnvkit_targets.bed",antitargets="cnvkit_antitargets.bed",genome=config['references'][pfamily]['CNVKITGENOME'],rname="pl:cnvkit_somatic"
    shell: "module load cnvkit/0.8; mkdir {output.dir}; cnvkit.py batch {input.tumor} -n {input.normal} --fasta {params.genome} --short-names --drop-low-coverage -g {params.access} --targets {params.targets} --antitargets {params.antitargets} --output-reference {output.dir}/{params.normalsample}.reference.cnn --output-dir {output.dir} --scatter; cnvkit.py call -o {output.calls} {output.dir}/{params.tumorsample}.cns; cnvkit.py gainloss -s {output.dir}/{params.tumorsample}.cns -t 0.3 --drop-low-coverage -o {output.gainloss} {output.dir}/{params.tumorsample}.cnr"