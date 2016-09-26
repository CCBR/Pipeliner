rule delly_somatic:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"
    output: deletion="delly_out/{x}_del.bcf",
            ins="delly_out/{x}_ins.bcf",
            dup="delly_out/{x}_dup.bcf",
            tra="delly_out/{x}_tra.bcf",
            inv="delly_out/{x}_inv.bcf"
    params: exclusions=config['references'][pfamily]['DELLY_EXCLUSIONS'],genome=config['references'][pfamily]['GENOME'],rname="pl:delly"
    shell: "module load delly; delly call -t DEL -x {params.exclusions} -o {output.deletion} -g {params.genome} {input.tumor} {input.normal}; delly call -t DUP -x {params.exclusions} -o {output.dup} -g {params.genome} {input.tumor} {input.normal}; delly call -t INV -x {params.exclusions} -o {output.inv} -g {params.genome} {input.tumor} {input.normal}; delly call -t TRA -x {params.exclusions} -o {output.tra} -g {params.genome} {input.tumor} {input.normal}; delly call -t INS -x {params.exclusions} -o {output.ins} -g {params.genome} {input.tumor} {input.normal}"