rule varscan:
       input:  lambda wildcards: config['project']['pairs'][wildcards.x][0]+".pileup.bam",
               lambda wildcards: config['project']['pairs'][wildcards.x][1]+".pileup.bam"
       output: snp="{x}.snp",
               indel= "{x}.indel",
               snpvcf="{x}.snp.vcf",
               indelvcf="{x}.indel.vcf"
       params: varscan=config['bin'][pfamily]['VARSCAN'],
               snpvcf="{x}.snp",indelvcf="{x}.indel",rname="pl:varscan"
       shell:  "{params.varscan} somatic {input[0]} {input[1]} --output-snp {output.snp} --output-indel {output.indel} --tumor-purity 0.85 --strand-filter 1;{params.varscan} somatic {input[0]} {input[1]} --output-snp {params.snpvcf} --output-indel {params.indelvcf} --output-vcf 1 --tumor-purity 0.85 --strand-filter 1"

