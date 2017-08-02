rule picard_rnaseqvar_markdups:
     input:  "{x}.sorted.bam"
     output: out = temp("{x}.dedup_unsplit.bam"),
             metrics = "{x}.sorted.txt"
     params: markdups=config['bin'][pfamily]['MARKDUPS'],rname="pl:rnamkdup"
     shell:  "{params.markdups} I={input} O={output.out} M={output.metrics} REMOVE_DUPLICATES=FALSE AS=TRUE PG='null' Validation_Stringency=LENIENT CREATE_INDEX=true"

