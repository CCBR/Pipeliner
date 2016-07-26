rule samtools_mpileup:
     input:   "{X}.recal.bam"
     output:  "{X}.pileup.bam"
     params:  genome=config['references']['GENOME'],,rname="pl:mpileup"
     shell: "samtools mpileup -q 1 -f {params.genome} {input} > {output}"

