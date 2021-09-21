rule novocraft_sort:
     input:  "{x}.fin.bam",
     output: temp("{x}.sorted.bam")
     params: novosort=config['bin'][pfamily]['NOVOSORT'],rname="pl:novosort"
     threads: 8
     shell:  "module load novocraft/3.08.02; novosort -t /lscratch/$SLURM_JOB_ID -m 100G --threads {threads} -s -i -o {output} {input};"


