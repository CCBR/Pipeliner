rule gatk_combine_somatic_gvcfs:
    input: "xa{batches}"
    output: temp("xa{batches}.gvcf")
    params: batch ="-l nodes=1:gpfs -q ccr",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],rname="pl:combgvcfs"
    run: 
       F=open(input[0],"r")
       fl=F.read().split()
       F.close()
       var=" --variant "+" --variant ".join(fl)
       cmd="module load GATK/3.8-0; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T CombineGVCFs -R {params.genome} -o "+output[0]+ var
       shell(cmd)