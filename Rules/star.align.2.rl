rule star_align_2:
     input:  file1="{x}.R1.trimmed.fastq.gz",file2="{x}.R2.trimmed.fastq.gz",length="QC/{x}_readlength.txt",tab=expand("{x}SJ.out.tab",x=samples)
     output: out1=temp("{x}.p2Aligned.sortedByCoord.out.bam"),out4="{x}.p2SJ.out.tab",out5="{x}.p2Log.final.out"
     params: rname='pl:star2p',prefix="{x}.p2",outsamunmapped=config['bin'][pfamily]['OUTSAMUNMAPPED'],wigtype=config['bin'][pfamily]['WIGTYPE'],wigstrand=config['bin'][pfamily]['WIGSTRAND'], gtffile=config['references'][pfamily]['FUSIONGTFFILE'], nbjuncs=config['bin'][pfamily]['NBJUNCS'],starref=config['references'][pfamily]['STARREF']
     threads: 32
     run:
        import glob
        rl=int(open(input.length).readlines()[0].strip())-1
        dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.starref+'*/',recursive=False))))
        bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
        cmd="module load STAR/2.5.2b; STAR --genomeDir {params.starref}"+str(bestdbrl)+" --readFilesIn {input.file1} {input.file2} --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {params.prefix} --outSAMunmapped {params.outsamunmapped} --sjdbFileChrStartEnd {input.tab} --sjdbGTFfile {params.gtffile} --outSAMtype BAM SortedByCoordinate"
        shell(cmd)