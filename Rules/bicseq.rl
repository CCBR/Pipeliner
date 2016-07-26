rule bicseq:
     input: bam="{x}.recal.bam"
     output: done="bicseq.{x}.done",dir="bicseq.{x}",description="sample_{x}"
     params: bicseq=config['bin'][pfamily]['BICSEQ'],genome=config['references'][pfamily]['GENOME'],rname="pl:bicseq"
     threads: 1
     shell: """
	samtools view -U BWA,{output.dir},N,N {input}; {params.bicseq} --bin_size=1000 --paired {configfile} {output.dir} {params.description};touch {output.done}
            """
