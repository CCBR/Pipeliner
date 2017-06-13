rule picard_rnaseqvarheaders:
     input:  "{x}.p2Aligned.sortedByCoord.out.bam",
     output: temp("{x}.sorted.bam")
     params: picard1=config['bin'][pfamily]['PICARD1'],rname="pl:rnahead"
     threads: 1
     run:  rgid=config['project']['units'][wildcards.x];rgpl=config['project']['rgid'][wildcards.x]['rgpl'];rgsm=config['project']['rgid'][wildcards.x]['rgsm'];rglb=config['project']['rgid'][wildcards.x]['rglb'];rgpu=config['project']['rgid'][wildcards.x]['rgpu'];rgcn=config['project']['rgid'][wildcards.x]['rgcn'];shell('{params.picard1} I={input} O={output} RGID={rgid} RGPL={rgpl} RGLB={rglb} RGPU={rgpu} RGSM={rgsm} RGCN={rgcn} RGDS={input} Validation_Stringency=LENIENT')
