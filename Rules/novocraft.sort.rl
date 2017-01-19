rule novocraft_sort:
     input:  "{x}.fin.bam",
     output: temp("{x}.sorted.bam")
     params: novosort=config['bin'][pfamily]['NOVOSORT'],rname="pl:novosort"
     threads: 2
     shell:  "module load novocraft/3.02.10;{params.novosort} -t /scratch --threads {threads} -s -i -o {output} {input};"


