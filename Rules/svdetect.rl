rule svdetect:
     input: lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam"
     output:"{x}.g.vcf"
     params: svdetect=config['bin'][pfamily]['SVDETECT'],genome=config['references'][pfamily]['GENOME'],rname="pl:svdetect"
     threads: 8
     shell: "{params.svdetect}"

