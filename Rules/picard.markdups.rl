rule picard_markdups:
     input:  "{x}.sorted.bam"
     output: out = temp("{x}.dedup.bam"),
             metrics = "{x}.sorted.txt"
     params: markdups=config['bin'][pfamily]['MARKDUPS'],rname="pl:mkdup"
     shell:  "{params.markdups} I={input} O={output.out} M={output.metrics} REMOVE_DUPLICATES=FALSE AS=TRUE PG='null' Validation_Stringency=LENIENT"

