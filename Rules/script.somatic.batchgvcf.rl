bch=str(BS) # combinegvcf batchsize
rule script_somatic_batchgvcf:
    input: expand(config['project']['workpath']+"/{s}.g.vcf",s=samples)
    output: expand(config['project']['workpath']+"/xa{x}",x=batches),"all.g.vcf"
    params: rname="pl:batchgvcf"
    run:
        # o=os.popen("rm all.g.vcf x*")
        #for s in samples:    
        #    o=os.popen("ls {0}.g.vcf >> all.g.vcf".format(s))
        #o=os.popen("split -l {} all.g.vcf".format(bch))
        o=os.popen("ls *.g.vcf>all.g.vcf;split -l {} all.g.vcf".format(bch))