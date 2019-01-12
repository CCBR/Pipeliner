rule make_target_files:
    input: expand("{s}.R1.trimmed.fastq.gz", s=samples)
    output: targets=config['project']['workpath']+"/exome_targets.bed",
             freectargets=config['project']['workpath']+"/freec_targets.bed",
    params: annot=config['project']['annotation'],bed=config['project']['targetspath'],access=config['references'][pfamily]['CNVKITACCESS'],gtf=config['references'][pfamily]['CNVKITANNO'],rname="pl:targets"
    shell: "perl Scripts/reformat_bed.pl {params.bed} {params.annot}"