rule pindel_germline:
    params: workdir=config['project']['workpath'],genome=config['references'][pfamily]['GENOME'],rname="pl:pindel"
    input: expand("{x}.recal.bam", x=samples)
    output: "pindel_out/pindel_calls",
    threads: 8
    shell: "mkdir -p pindel_out; perl Scripts/make_pindel_config.pl; module load pindel/0.2.5b8; pindel -i pindel_out/pindel_config -f (params.genome) -o {output} --number_of_threads {threads}"