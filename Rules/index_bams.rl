rule index_bams:
     input:  "{x}.recal.bam"
     output: out1="{x}.recal.bam.bai", out2="{x}.recal.bai"
     params: markdups=config['bin'][pfamily]['INDEXBAM'],rname="pl:indexbam"
     shell:  "module load samtools/1.9; samtools index {input}; cp {output.out1} {output.out2}"

