rule svaba_wgs_somatic:
     input:  normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
             tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
     output: "svaba_out/{x}.log"
     params: tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],genome=config['references'][pfamily]['GENOME'],rname="svaba"
     threads: 32
     shell: "mkdir -p svaba_out; cd svaba_out; module load gcc; /data/CCBR/apps/svaba/bin/svaba run -t ../{input.tumor} -n ../{input.normal} -a {params.normalsample}+{params.tumorsample} -p {threads} -G {params.genome}"