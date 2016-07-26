rule picard_markdups:
     input:  "{x}.sorted.bam"
     output: out = "{x}.dedup.bam",
             metrics = "{x}.sorted.txt"
     params: markdups=config['bin']['MARKDUPS'],rname="pl:mkdup"
     shell:  "{params.markdups} I={input} O={output.out} M={output.metrics} REMOVE_DUPLICATES=TRUE AS=TRUE PG='null' Validation_Stringency=LENIENT"

