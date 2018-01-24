rule manta_somatic:
     input:  normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
             tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
     output: vcf="manta_out/{x}/results/variants/candidateSV.vcf.gz",
     params: tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0], gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:manta_somatic"
     threads: 8
     shell: "mkdir -p manta_out/{params.normalsample}+{params.tumorsample}; module load manta/1.2.0; module load python/2.7; configManta.py --bam={input.normal} --tumorBam={input.tumor} --referenceFasta {params.genome} --exome --runDir manta_out/{params.normalsample}+{params.tumorsample}; manta_out/{params.normalsample}+{params.tumorsample}/runWorkflow.py -m local -j {threads} -g 12"