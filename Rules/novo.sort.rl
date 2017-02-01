rule novosort_sort:
     input:  "{x}.fin.bam"
     output: temp("{x}.sorted.bam")
     params: novosort=config['bin'][pfamily]['NOVOSORT'],rname="pl:novosort"
     threads: 8
     shell:  "module load novocraft;{params.novosort} --threads {threads} -t /scratch -m 100G -s -i -o {output} {input};"


