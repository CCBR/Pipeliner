rule delly_germline:
    params: exclusions=config['references'][pfamily]['DELLY_EXCLUSIONS'],workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pl:delly"
    input: expand("{x}.recal.bam", x=samples)
    output: deletion=config['project']['workpath']+"/delly_out/deletions.bcf",
            ins=config['project']['workpath']+"/delly_out/insertions.bcf",
            dup=config['project']['workpath']+"/delly_out/duplications.bcf",
            tra=config['project']['workpath']+"/delly_out/translocations.bcf",
            inv=config['project']['workpath']+"/delly_out/inversions.bcf"
    shell: "mkdir -p delly_out; module load delly/0.7.6; delly call -t DEL -x {params.exclusions} -o {output.deletion} -g {params.genome} {input}; delly call -t DUP -x {params.exclusions} -o {output.dup} -g {params.genome} {input}; delly call -t INV -x {params.exclusions} -o {output.inv} -g {params.genome} {input}; delly call -t TRA -x {params.exclusions} -o {output.tra} -g {params.genome} {input}; delly call -t INS -x {params.exclusions} -o {output.ins} -g {params.genome} {input}"