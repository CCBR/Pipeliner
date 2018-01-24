rule manta_somatic_tumoronly:
     input:  tumor=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam"
     output: vcf="manta_out/{x}/results/variants/candidateSV.vcf.gz",
     params: tumorsample=lambda wildcards: config['project']['units'][wildcards.x],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:manta"
     threads: 8
     shell: "mkdir -p manta_out/{params.tumorsample}; module load manta/1.2.0; module load python/2.7; configManta.py --tumorBam={input.tumor} --referenceFasta {params.genome} --exome --runDir manta_out/{params.tumorsample}; manta_out/{params.tumorsample}/runWorkflow.py -m local -j {threads} -g 12"