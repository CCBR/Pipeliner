rule delly_germline_del:
    params: exclusions=config['references'][pfamily]['DELLY_EXCLUSIONS'],workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pl:delly"
    input: expand("{x}.recal.bam", x=samples)
    output: deletion=config['project']['workpath']+"/delly_out/deletions.bcf",
    shell: "mkdir -p delly_out; module load delly/0.7.6; delly call -t DEL -x {params.exclusions} -o {output.deletion} -g {params.genome} {input}; module load samtools; bcftools convert -O v -o deletions.vcf {output.deletion}"

rule delly_germline_ins:
    params: exclusions=config['references'][pfamily]['DELLY_EXCLUSIONS'],workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pl:delly"
    input: expand("{x}.recal.bam", x=samples)
    output: ins=config['project']['workpath']+"/delly_out/insertions.bcf",
    shell: "mkdir -p delly_out; module load delly/0.7.6; delly call -t INS -x {params.exclusions} -o {output.ins} -g {params.genome} {input}; module load samtools; bcftools convert -O v -o insertions.vcf {output.ins}"

rule delly_germline_dup:
    params: exclusions=config['references'][pfamily]['DELLY_EXCLUSIONS'],workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pl:delly"
    input: expand("{x}.recal.bam", x=samples)
    output: dup=config['project']['workpath']+"/delly_out/duplications.bcf",
    shell: "mkdir -p delly_out; module load delly/0.7.6; delly call -t DUP -x {params.exclusions} -o {output.dup} -g {params.genome} {input}; module load samtools; bcftools convert -O v -o duplications.vcf {output.dup}"

rule delly_germline_trans:
    params: exclusions=config['references'][pfamily]['DELLY_EXCLUSIONS'],workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pl:delly"
    input: expand("{x}.recal.bam", x=samples)
    output: tra=config['project']['workpath']+"/delly_out/translocations.bcf",
    shell: "mkdir -p delly_out; module load delly/0.7.6; delly call -t TRA -x {params.exclusions} -o {output.tra} -g {params.genome} {input}; module load samtools; bcftools convert -O v -o translocations.vcf {output.tra}"

rule delly_germline_inv:
    params: exclusions=config['references'][pfamily]['DELLY_EXCLUSIONS'],workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pl:delly"
    input: expand("{x}.recal.bam", x=samples)
    output: inv=config['project']['workpath']+"/delly_out/inversions.bcf"
    shell: "mkdir -p delly_out; module load delly/0.7.6; delly call -t INV -x {params.exclusions} -o {output.inv} -g {params.genome} {input}; module load samtools; bcftools convert -O v -o inversions.vcf {output.inv}"