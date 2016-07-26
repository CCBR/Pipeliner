rule star_align_2:
     input:  "{x}.R1."+config['project']['filetype'],
             "{x}.R2."+config['project']['filetype']
     output: temp("{x}.p2.bam")
     params: genomeDir=config['references']['genomeDir2'],
             gtf=config['references']['gencodeGtf'],
             adapter1=config['references']['ADAPTER1'],
             adapter2=config['references']['ADAPTER2'],rname="pl:star2"
     threads: 8
     version: "1.0"
     shell: "module load STAR;  STAR --genomeDir {params.genomeDir2} --readFilesIn {input} --readFilesCommand zcat --runThreadN {threads} --sjdbGTFfile {params.genomeDir}/{params.gtf} --sjdbOverhang 99 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --clip3pAdapterSeq {params.adapter1} {params.adapter2} --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate"
     
#      --outFileNamePrefix $3/pass2/

