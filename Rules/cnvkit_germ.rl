rule cnvkit_germ:
    input: expand("{x}.recal.bam", x=samples)
    output: heatmap="cnvkit_out/cnvkit_heatmap.pdf"
    threads: 12
    params: cnvkitgenome=config['references'][pfamily]['CNVKITREF'],rname="pl:cnvkit"
    shell: "mkdir -p cnvkit_out; cd cnvkit_out; module load cnvkit/0.9.1; cnvkit.py batch --scatter --method wgs -r {params.cnvkitgenome} -p {threads} ../*.recal.bam; cnvkit.py heatmap -d -o cnvkit_heatmap.pdf *.cns"