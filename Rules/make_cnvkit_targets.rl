rule make_target_files:
    input: expand("{s}.R1.trimmed.fastq.gz", s=samples)
    output: targets=config['project']['workpath']+"/exome_targets.bed",
            cnvkittargets=config['project']['workpath']+"/cnvkit_targets.bed",
            cnvkitantitargets=config['project']['workpath']+"/cnvkit_antitargets.bed"
    params: access=config['references'][pfamily]['CNVKITACCESS'],gtf=snpsites=config['references'][pfamily]['GTFFILE'],rname="pl:targets"
    shell: "module load cnvkit/0.8.5; cnvkit.py target {params.access} --split --short-names --annotate {params.gtf} -o targets.bed -a 1000"