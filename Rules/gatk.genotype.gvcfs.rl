rule gatk_genotype_gvcfs:
    input: expand("xa{x}.gvcf",x=batches)
    output: "combined.vcf"
    params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
    threads: 12
    run:
      fl=os.popen("ls x*.gvcf").read().split()      
      var=" --variant "+" --variant ".join(fl)
      cmd="module load GATK/3.6; GATK -m 120G GenotypeGVCFs -R {params.genome} --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -nt {threads}"+var
      print(cmd)
      shell(cmd)

