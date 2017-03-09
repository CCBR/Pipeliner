rule picard_headers:
     input:  "{x}.bam",
     output: temp("{x}.fin.bam")
     params: picard1=config['bin'][pfamily]['PICARD1'],rname="pl:headers"
     run:  rgid=config['project']['units'][wildcards.x];rgpl=config['project']['platform'];shell('{params.picard1} I={input} O={output} RGID={rgid} RGPL={rgpl} RGLB={rgid} RGPU=1_1 RGSM={rgid} RGCN=CCBR RGDS={input} Validation_Stringency=LENIENT')
