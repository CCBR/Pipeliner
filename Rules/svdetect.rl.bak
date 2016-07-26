rule svdetect:
     input: lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam"
     output:"{x}.g.vcf"
     params: svdetect=config['bin']['SVDETECT'],genome=config['references']['GENOME'],rname="pl:svdetect"
     threads: 8
     shell: "{params.svdetect}"

