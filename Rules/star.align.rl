rule star_align_1:
     input:  "{x}.R1."+config['project']['filetype'],
             "{x}.R2."+config['project']['filetype']
     output: temp("{x}.p1.bam")
     params: genomeDir=config['references'][pfamily]['genomeDir'],
             gtf=config['references'][pfamily]['gencodeGtf'],
             adapter1=config['references'][pfamily]['ADAPTER1'],
             adapter2=config['references'][pfamily]['ADAPTER2'],rname="pl:star"
     threads: 8
     version: "1.0"
     shell: "module load STAR; STAR --genomeDir {params.genomeDir} --readFilesIn {input} --readFilesCommand zcat --runThreadN {threads} --sjdbGTFfile {params.genomeDir}/{params.gencodeGtf} --sjdbOverhang 99 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterType BySJout --clip3pAdapterSeq {params.adapter1} {params.adapter2} --outSAMtype BAM SortedByCoordinate --chimSegmentMin 25 --outStd SAM > {output}"
# --outFileNamePrefix $3 
