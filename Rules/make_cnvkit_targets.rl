rule make_cnvkit_targets:
    input: expand("{s}.R1.trimmed.fastq.gz", s=samples)
    output: "targets.bed",
    params: access=config['references'][pfamily]['CNVKITACCESS'],gtf=config['references'][pfamily]['CNVKITANNO'],rname="pl:targets"
    shell: "module load cnvkit/0.8.5; cnvkit.py target {params.access} --split --short-names --annotate {params.gtf} -o {output} -a 1000"