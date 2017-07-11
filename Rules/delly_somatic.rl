rule delly_somatic_del:
    params: exclusions=config['references'][pfamily]['DELLY_EXCLUSIONS'],workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pl:dellydel"
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"
    output: deletion=config['project']['workpath']+"/delly_out/{x}_del.bcf",delvcf=config['project']['workpath']+"/delly_out/{x}_del.vcf"
    shell: "module load delly; delly call -t DEL -x {params.exclusions} -o {output.deletion} -g {params.genome} {input.tumor} {input.normal}; module load samtools; bcftools convert -O v -o {output.delvcf} {output.deletion}"

rule delly_somatic_ins:
    params: exclusions=config['references'][pfamily]['DELLY_EXCLUSIONS'],workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pl:dellyins"
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"
    output: ins=config['project']['workpath']+"/delly_out/{x}_ins.bcf",insvcf=config['project']['workpath']+"/delly_out/{x}_ins.vcf"
    shell: "module load delly; delly call -t INS -x {params.exclusions} -o {output.ins} -g {params.genome} {input.tumor} {input.normal}; module load samtools; bcftools convert -O v -o {output.insvcf} {output.ins}"

rule delly_somatic_dup:
    params: exclusions=config['references'][pfamily]['DELLY_EXCLUSIONS'],workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pl:dellydup"
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"
    output: dup=config['project']['workpath']+"/delly_out/{x}_dup.bcf",dupvcf=config['project']['workpath']+"/delly_out/{x}_dup.vcf"
    shell: "module load delly; delly call -t DUP -x {params.exclusions} -o {output.dup} -g {params.genome} {input.tumor} {input.normal}; module load samtools; bcftools convert -O v -o {output.dupvcf} {output.dup}"

rule delly_somatic_trans:
    params: exclusions=config['references'][pfamily]['DELLY_EXCLUSIONS'],workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pl:dellytrans"
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"
    output: tra=config['project']['workpath']+"/delly_out/{x}_tra.bcf",transvcf=config['project']['workpath']+"/delly_out/{x}_trans.vcf"
    shell: "module load delly; delly call -t TRA -x {params.exclusions} -o {output.tra} -g {params.genome} {input.tumor} {input.normal}; module load samtools; bcftools convert -O v -o {output.transvcf} {output.tra}"

rule delly_somatic_inv:
    params: exclusions=config['references'][pfamily]['DELLY_EXCLUSIONS'],workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pl:dellyinv"
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"
    output: inv=config['project']['workpath']+"/delly_out/{x}_inv.bcf",invvcf=config['project']['workpath']+"/delly_out/{x}_inv.vcf"
    shell: "module load delly; delly call -t INV -x {params.exclusions} -o {output.inv} -g {params.genome} {input.tumor} {input.normal}; module load samtools; bcftools convert -O v -o {output.invvcf} {output.inv}"