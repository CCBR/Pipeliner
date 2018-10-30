rule svaba_wgs_somatic:
     input:  normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
             tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
             normalbai=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam.bai",
             tumorbai=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam.bai"
     output: "svaba_out/{x}.log"
     params: tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],genome=config['references'][pfamily]['CNVKITGENOME'],rname="svaba"
     threads: 32
     shell: "mkdir -p svaba_out; cd svaba_out; module load gcc; /data/CCBR_Pipeliner/db/PipeDB/bin/svaba/svaba/bin/svaba run -t ../{input.tumor} -n ../{input.normal} -a {params.normalsample}+{params.tumorsample} -p {threads} -G {params.genome}"