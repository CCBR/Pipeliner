rule gatk_combine_gvcfs:
    input: "xa{batches}"
    output:"xa{batches}.gvcf"
    params: batch ="-l nodes=1:gpfs -q ccr",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],rname="pl:combgvcfs"
    run: 
       F=open(input[0],"r")
       fl=F.read().split()
       F.close()
       var=" --variant "+" --variant ".join(fl)
       cmd="{params.gatk} -R {params.genome} -T CombineGVCFs -o "+output[0]+ var
       shell(cmd)


