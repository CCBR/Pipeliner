rule cnvnator_wgs:
     input: lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam"
     output:"{x}.g.vcf"
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="cnvnator"
     shell: "module load cnvnator/0.3.3; "


