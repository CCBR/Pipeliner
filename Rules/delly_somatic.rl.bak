rule delly_somatic:
    input: lambda wildcards: list(config['project']['pairs'].values())[wildcards.group],
    output: del="delly/{group[1]}.del.bcf",
            ins="delly/{group[1]}.ins.bcf",
            dup="delly/{group[1]}.dup.bcf",
            tra="delly/{group[1]}.tra.bcf",
            inv="delly/{group[1]}.inv.bcf"
    params: exclusions=config['references']['DELLY_EXCLUSIONS'],genome=config['references']['GENOME'],rname="pl:delly"
    shell: "module load delly; delly call -t DEL -x {params.exclusions} -o {output.del} -g {params.genome} {group[1]}.realign.bam {group[0]}.realign.bam; delly call -t DUP -x {params.exclusions} -o {output.dup} -g {params.genome} {group[1]}.realign.bam {group[0]}.realign.bam; delly call -t INV -x {params.exclusions} -o {output.inv} -g {params.genome} {group[1]}.realign.bam {group[0]}.realign.bam; delly call -t TRA -x {params.exclusions} -o {output.tra} -g {params.genome} {group[1]}.realign.bam {group[0]}.realign.bam; delly call -t INS -x {params.exclusions} -o {output.ins} -g {params.genome} {group[1]}.realign.bam {group[0]}.realign.bam"
