rule samtools_sort:
     input:  "{x}.p2Aligned.out.bam",
     output: temp("{x}.p2Aligned.sortedByCoord.out.bam")
     params: novosort=config['bin'][pfamily]['NOVOSORT'],rname="pl:sort"
     threads: 8
     shell:  "module load samtools/1.9; samtools sort -@ {threads} -o {output} {input};"


