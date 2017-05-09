rule cnvkit_germ:
    input: expand("{x}.recal.bam", x=samples)
    output: heatmap="cnvkit_out/germline_cnvkit.heatmap"
    threads: 12
    params: cnvkitgenome=config['references'][pfamily]['CNVKITGENOME'],rname="pl:cnvkit"
    shell: "mkdir -p cnvkit_out; cd cnvkit_out; module load cnvkit/0.8.5; cnvkit.py batch --scatter --method wgs -r {params.cnvkitgenome} -p {threads} ../*.recal.bam; cnvkit.py heatmap -d -o germline_cnvkit.heatmap *.cns"