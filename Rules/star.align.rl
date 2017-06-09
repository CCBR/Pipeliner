rule star_align_1:
     input:  "{x}.R1.trimmed.fastq.gz","{x}.R2.trimmed.fastq.gz"
     output: out1=temp("{x}.p2.Aligned.sortedByCoord.out.bam"),out2="{x}.p2.ReadsPerGene.out.tab",out3="{x}.p2.Aligned.toTranscriptome.out.bam",out4="{x}.p2.SJ.out.tab",out5="{x}.p2.Log.final.out"
     params: rname='pl:star1p',prefix="{x}",starref=config['references'][pfamily]['STARREF']+config['project']["SJDBOVERHANG"]
     threads: 32
     shell: "module load STAR/2.5.2a; STAR --genomeDir {params.starref} --readFilesIn {input.file1} {input.file2} --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {params.prefix}. --outSAMtype BAM Unsorted"