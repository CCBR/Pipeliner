rule script_bzip2:
    input:  "{x}.fastq"
    output: "{x}.fastq.bz2"
    params: rname="pl:bzip2"
    shell:  "bzip2 {input} > {output}"

