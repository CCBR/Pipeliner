rule star_align_1:
     input:  file1="{x}.R1.trimmed.fastq.gz",file2="{x}.R2.trimmed.fastq.gz",file3="QC/{x}.R1.trimmed_fastqc.html",file4="QC/{x}.R2.trimmed_fastqc.html",length="QC/{x}_readlength.txt"
     output: out1=temp("{x}Aligned.out.bam"),out4="{x}SJ.out.tab",out5="{x}Log.progress.out"
     params: rname='pl:star1p',prefix="{x}",starref=config['references'][pfamily]['STARREF']+config['project']["SJDBOVERHANG"]
     threads: 32
     run:
        import glob,json
        rl=int(open("{input.length}").readlines()[0].strip())-1
        dbrl=sorted(list(map(lambda x:int(re.findall("genes-(\d+)",x)[0]),glob.glob(params.stardir+'*/',recursive=False))))
        bestdbrl=next(x[1] for x in enumerate(dbrl) if x[1] >= rl)
        cmd="module load STAR/2.5.2b; STAR --genomeDir {params.starref}"+str(bestdbrl)+" --readFilesIn {input.file1} {input.file2} --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {params.prefix} --outSAMtype BAM Unsorted"
        shell(cmd)
        # add sjdboverhand and stardir back to run.json
        A=open("run.json",'r')
        a=eval(A.read())
        A.close()
        config=dict(a.items())
        config['project']['SJDBOVERHANG']=str(bestdbrl)
        config['project']["STARDIR"]= config['references'][pfamily]['STARDIR']+str(bestdbrl)
        with open('run.json','w') as F:
          json.dump(config, F, sort_keys = True, indent = 4,ensure_ascii=False)
        F.close()