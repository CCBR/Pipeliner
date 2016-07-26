rule novosort_sort:
     input:  "{x}.fin.bam"
     output: temp("{x}.sorted.bam")
     params: novosort=config['bin'][pfamily]['NOVOSORT'],rname="pl:novosort"
     shell:  "module load novocraft;{params.novosort} -s -i -o {output} {input};"


