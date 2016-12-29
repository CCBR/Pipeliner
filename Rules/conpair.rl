rule conpair:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           dir="conpair_out"
    output: estimate="conpair_out/{x}.conpair",
    params: tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],genome=config['references'][pfamily]['GENOME'],markers=config['references'][pfamily]['CONPAIRMARKERS'],rname="pl:conpair"
    shell: "module load conpair; run_gatk_pileup_for_sample.py -B {input.normal} -R {params.genome} -O conpair_out/{params.normalsample}.pileup; run_gatk_pileup_for_sample.py -B {input.tumor} -R {params.genome} -O conpair_out/{params.tumorsample}.pileup; estimate_tumor_normal_contamination.py -T conpair_out/{params.tumorsample}.pileup -N conpair_out/{params.normalsample}.pileup -O {output} -M {params.markers}"