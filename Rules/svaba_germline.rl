rule svaba_germline:
     input: expand("{x}.recal.bam", x=samples)
     output: config['project']['workpath']+"/svaba_out/svaba.log"
     params: genome=config['references'][pfamily]['CNVKITGENOME'],rname="svaba"
     threads: 32
     run:
        fl=os.popen("ls *.recal.bam").read().split()      
        var=" -t ../"+" -t ../".join(fl)
        cmd="mkdir -p svaba_out; cd svaba_out; module load gcc; /data/CCBR_Pipeliner/db/PipeDB/bin/svaba/svaba/bin/svaba run -a svaba -p {threads} --germline -G {params.genome}"+var
        print(cmd)
        shell(cmd)