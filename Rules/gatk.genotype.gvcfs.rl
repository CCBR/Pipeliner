if config['project']['annotation'] == "hg19":

  rule gatk_genotype_gvcfs_1:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_1.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 1 -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_2:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_2.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 2 -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_3:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_3.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 3 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_4:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_4.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 4 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_5:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_5.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 5 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_6:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_6.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 6 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_7:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_7.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 7 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_8:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_8.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 8 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_9:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_9.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 9 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_10:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_10.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 10 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_11:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_11.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 11 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_12:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_12.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 12 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_13:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_13.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 13 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_14:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_14.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 14 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_15:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_15.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 15 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_16:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_16.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 16 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_17:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_17.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 17 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_18:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_18.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 18 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_19:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_19.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 19 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_20:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_20.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 20 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_21:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_21.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 21 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_22:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_22.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 22 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_X:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_X.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L X -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_Y:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_Y.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L Y -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_MT:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_MT.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L MT -nt {threads}"+var
        print(cmd)
        shell(cmd)


elif config['project']['annotation'] == "hg38":

  rule gatk_genotype_gvcfs_1:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_1.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr1 -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_2:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_2.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr2 -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_3:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_3.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr3 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_4:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_4.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr4 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_5:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_5.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr5 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_6:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_6.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr6 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_7:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_7.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr7 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_8:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_8.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr8 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_9:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_9.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr9 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_10:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_10.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr10 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_11:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_11.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr11 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_12:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_12.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr12 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_13:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_13.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr13 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_14:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_14.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr14 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_15:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_15.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr15 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_16:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_16.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr16 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_17:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_17.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr17 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_18:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_18.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr18 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_19:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_19.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr19 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_20:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_20.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr20 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_21:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_21.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr21 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_22:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_22.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr22 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_X:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_X.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chrX -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_Y:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_Y.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chrY -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_MT:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_MT.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chrM -nt {threads}"+var
        print(cmd)
        shell(cmd)


elif config['project']['annotation'] == "mm10":

  rule gatk_genotype_gvcfs_1:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_1.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr1 -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_2:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_2.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr2 -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_3:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_3.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr3 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_4:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_4.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr4 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_5:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_5.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr5 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_6:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_6.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr6 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_7:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_7.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr7 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_8:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_8.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr8 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_9:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_9.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr9 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_10:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_10.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr10 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_11:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_11.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr11 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_12:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_12.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr12 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_13:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_13.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr13 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_14:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_14.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr14 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_15:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_15.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr15 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_16:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_16.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr16 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_17:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_17.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr17 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_18:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_18.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr18 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_19:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_19.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr19 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_X:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_X.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chrX -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_Y:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_Y.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chrY -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_MT:
      input: expand("{x}.g.vcf",x=samples)
      output: temp("combined_MT.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-0; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chrM -nt {threads}"+var
        print(cmd)
        shell(cmd)