rule make_target_files:
    input: config['project'][targetspath]
    output: targets="exome_targets.bed",
            cnvkittargets="cnvkit_targets.bed",
            cnvkitantitargets="cnvkit_antitargets.bed"
    params: bed=config['project'][targetspath],rname="pl:targets"
    shell: "perl Scripts/reformat_bed.pl {input}; module load cnvkit; cnvkit.py target --split -o {output.cnvkittargets} {input}; cnvkit.py antitarget -o {output.cnvkitantitargets} {input}"