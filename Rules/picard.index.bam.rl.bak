rule picard_index_bam:
     input:  "{x}.bam"
     output: out = "{x}.bai"
     params: markdups=config['bin']['INDEXBAM'],rname="pl:indexbam"
     shell:  "{params.markdups} I={input} O={output.out}"

