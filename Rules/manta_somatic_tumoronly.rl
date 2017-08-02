rule manta_somatic_tumoronly:
     input:  lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam"
     output: vcf="manta_out/{x}/results/variants/candidateSV.vcf.gz",
     params: tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:manta"
     threads: 8
     shell: "mkdir -p manta_out/{params.tumorsample}; module load manta/1.1.0; module load python/2.7; configManta.py --tumorBam={input.tumor} --referenceFasta {params.genome} --exome --runDir manta_out/{params.normalsample}+{params.tumorsample}; manta_out/{params.normalsample}+{params.tumorsample}/runWorkflow.py -m local -j {threads} -g 12"