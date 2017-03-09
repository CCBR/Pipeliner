rule manta_somatic:
    input:   normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
             tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
             dir=config['project']['workpath']+"/manta_out"
     output: vcf=config['project']['workpath']+"/manta_out/{x}_candidateSV.vcf.gz",
             outdir=config['project']['workpath']+"/manta_out/{x}"
     params: tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0], gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:manta_somatic"
     threads: 8
     shell: "mkdir -p {output.outdir}; module load manta; module load python/2.7; configManta.py --bam={input.normal} --tumorBam={input.tumor} --referenceFasta {params.genome} --exome --runDir {output.outdir}; {output.outdir}/runWorkflow.py -m local -j {threads} -g 12; mv {input.dir}/{params.normalsample}+{params.tumorsample}/results/variants/candidateSV.vcf.gz {output.vcf}"