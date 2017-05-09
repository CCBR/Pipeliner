rule cnvkit_germ:
    input: bams=expand(lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam")
    output: cnr=expand("cnvkit_out/"+lambda wildcards: config['project']['units'][wildcards.x]+".cnr"),
            cns=expand("cnvkit_out/"+lambda wildcards: config['project']['units'][wildcards.x]+".cns"),
            heatmap="cnvkit/germline_cnvkit.heatmap"
    threads: 12
    params: cnvkitgenome=config['references'][pfamily]['CNVKITGENOME'],rname="pl:cnvkit"
    shell: "mkdir -p cnvkit_out; cd cnvkit_out; module load cnvkit/0.8.5; cnvkit.py batch --scatter --method wgs -r {params.cnvkitgenome} -p {threads} ../*.recal.bam; cnvkit.py heatmap -d -o cnvkit_heatmap.pdf *.cns"