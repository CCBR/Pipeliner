if config['project']['annotation'] == "hg19":

  rule gatk_genotype_gvcfs_1:
      input: expand("{x}_1.g.vcf",x=samples)
      output: temp("combined_1.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_1.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 1 -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_2:
      input: expand("{x}_2.g.vcf",x=samples)
      output: temp("combined_2.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_2.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 2 -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_3:
      input: expand("{x}_3.g.vcf",x=samples)
      output: temp("combined_3.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_3.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 3 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_4:
      input: expand("{x}_4.g.vcf",x=samples)
      output: temp("combined_4.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_4.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 4 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_5:
      input: expand("{x}_5.g.vcf",x=samples)
      output: temp("combined_5.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_5.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 5 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_6:
      input: expand("{x}_6.g.vcf",x=samples)
      output: temp("combined_6.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_6.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 6 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_7:
      input: expand("{x}_7.g.vcf",x=samples)
      output: temp("combined_7.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_7.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 7 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_8:
      input: expand("{x}_8.g.vcf",x=samples)
      output: temp("combined_8.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_8.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 8 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_9:
      input: expand("{x}_9.g.vcf",x=samples)
      output: temp("combined_9.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_9.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 9 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_10:
      input: expand("{x}_10.g.vcf",x=samples)
      output: temp("combined_10.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_10.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 10 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_11:
      input: expand("{x}_11.g.vcf",x=samples)
      output: temp("combined_11.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_11.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 11 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_12:
      input: expand("{x}_12.g.vcf",x=samples)
      output: temp("combined_12.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_12.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 12 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_13:
      input: expand("{x}_13.g.vcf",x=samples)
      output: temp("combined_13.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_13.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 13 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_14:
      input: expand("{x}_14.g.vcf",x=samples)
      output: temp("combined_14.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_14.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 14 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_15:
      input: expand("{x}_15.g.vcf",x=samples)
      output: temp("combined_15.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_15.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 15 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_16:
      input: expand("{x}_16.g.vcf",x=samples)
      output: temp("combined_16.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_16.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 16 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_17:
      input: expand("{x}_17.g.vcf",x=samples)
      output: temp("combined_17.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_17.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 17 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_18:
      input: expand("{x}_18.g.vcf",x=samples)
      output: temp("combined_18.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_18.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 18 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_19:
      input: expand("{x}_19.g.vcf",x=samples)
      output: temp("combined_19.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_19.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 19 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_20:
      input: expand("{x}_20.g.vcf",x=samples)
      output: temp("combined_20.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_20.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 20 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_21:
      input: expand("{x}_21.g.vcf",x=samples)
      output: temp("combined_21.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_21.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 21 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_22:
      input: expand("{x}_22.g.vcf",x=samples)
      output: temp("combined_22.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_22.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L 22 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_X:
      input: expand("{x}_X.g.vcf",x=samples)
      output: temp("combined_X.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_X.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L X -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_Y:
      input: expand("{x}_Y.g.vcf",x=samples)
      output: temp("combined_Y.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_Y.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L Y -nt {threads}"+var
        print(cmd)
        shell(cmd)

elif config['project']['annotation'] == "hg38":

  rule gatk_genotype_gvcfs_1:
      input: expand("{x}_1.g.vcf",x=samples)
      output: temp("combined_1.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_1.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr1 -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_2:
      input: expand("{x}_2.g.vcf",x=samples)
      output: temp("combined_2.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_2.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr2 -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_3:
      input: expand("{x}_3.g.vcf",x=samples)
      output: temp("combined_3.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_3.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr3 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_4:
      input: expand("{x}_4.g.vcf",x=samples)
      output: temp("combined_4.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_4.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr4 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_5:
      input: expand("{x}_5.g.vcf",x=samples)
      output: temp("combined_5.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_5.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr5 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_6:
      input: expand("{x}_6.g.vcf",x=samples)
      output: temp("combined_6.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_6.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr6 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_7:
      input: expand("{x}_7.g.vcf",x=samples)
      output: temp("combined_7.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_7.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr7 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_8:
      input: expand("{x}_8.g.vcf",x=samples)
      output: temp("combined_8.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_8.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr8 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_9:
      input: expand("{x}_9.g.vcf",x=samples)
      output: temp("combined_9.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_9.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr9 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_10:
      input: expand("{x}_10.g.vcf",x=samples)
      output: temp("combined_10.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_10.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr10 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_11:
      input: expand("{x}_11.g.vcf",x=samples)
      output: temp("combined_11.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_11.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr11 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_12:
      input: expand("{x}_12.g.vcf",x=samples)
      output: temp("combined_12.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_12.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr12 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_13:
      input: expand("{x}_13.g.vcf",x=samples)
      output: temp("combined_13.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_13.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr13 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_14:
      input: expand("{x}_14.g.vcf",x=samples)
      output: temp("combined_14.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_14.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr14 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_15:
      input: expand("{x}_15.g.vcf",x=samples)
      output: temp("combined_15.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_15.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr15 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_16:
      input: expand("{x}_16.g.vcf",x=samples)
      output: temp("combined_16.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_16.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr16 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_17:
      input: expand("{x}_17.g.vcf",x=samples)
      output: temp("combined_17.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_17.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr17 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_18:
      input: expand("{x}_18.g.vcf",x=samples)
      output: temp("combined_18.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_18.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr18 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_19:
      input: expand("{x}_19.g.vcf",x=samples)
      output: temp("combined_19.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_19.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr19 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_20:
      input: expand("{x}_20.g.vcf",x=samples)
      output: temp("combined_20.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_20.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr20 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_21:
      input: expand("{x}_21.g.vcf",x=samples)
      output: temp("combined_21.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_21.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr21 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_22:
      input: expand("{x}_22.g.vcf",x=samples)
      output: temp("combined_22.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_22.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr22 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_X:
      input: expand("{x}_X.g.vcf",x=samples)
      output: temp("combined_X.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_X.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chrX -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_Y:
      input: expand("{x}_Y.g.vcf",x=samples)
      output: temp("combined_Y.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_Y.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chrY -nt {threads}"+var
        print(cmd)
        shell(cmd)


elif config['project']['annotation'] == "mm10":

  rule gatk_genotype_gvcfs_1:
      input: expand("{x}_1.g.vcf",x=samples)
      output: temp("combined_1.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_1.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr1 -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_2:
      input: expand("{x}_2.g.vcf",x=samples)
      output: temp("combined_2.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_2.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr2 -nt {threads}"+var
        print(cmd)
        shell(cmd)

  rule gatk_genotype_gvcfs_3:
      input: expand("{x}_3.g.vcf",x=samples)
      output: temp("combined_3.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_3.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr3 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_4:
      input: expand("{x}_4.g.vcf",x=samples)
      output: temp("combined_4.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_4.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr4 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_5:
      input: expand("{x}_5.g.vcf",x=samples)
      output: temp("combined_5.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_5.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr5 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_6:
      input: expand("{x}_6.g.vcf",x=samples)
      output: temp("combined_6.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_6.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr6 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_7:
      input: expand("{x}_7.g.vcf",x=samples)
      output: temp("combined_7.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_7.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr7 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_8:
      input: expand("{x}_8.g.vcf",x=samples)
      output: temp("combined_8.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_8.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr8 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_9:
      input: expand("{x}_9.g.vcf",x=samples)
      output: temp("combined_9.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_9.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr9 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_10:
      input: expand("{x}_10.g.vcf",x=samples)
      output: temp("combined_10.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_10.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr10 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_11:
      input: expand("{x}_11.g.vcf",x=samples)
      output: temp("combined_11.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_11.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr11 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_12:
      input: expand("{x}_12.g.vcf",x=samples)
      output: temp("combined_12.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_12.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr12 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_13:
      input: expand("{x}_13.g.vcf",x=samples)
      output: temp("combined_13.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_13.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr13 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_14:
      input: expand("{x}_14.g.vcf",x=samples)
      output: temp("combined_14.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_14.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr14 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_15:
      input: expand("{x}_15.g.vcf",x=samples)
      output: temp("combined_15.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_15.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr15 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_16:
      input: expand("{x}_16.g.vcf",x=samples)
      output: temp("combined_16.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_16.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr16 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_17:
      input: expand("{x}_17.g.vcf",x=samples)
      output: temp("combined_17.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_17.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr17 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_18:
      input: expand("{x}_18.g.vcf",x=samples)
      output: temp("combined_18.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_18.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr18 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_19:
      input: expand("{x}_19.g.vcf",x=samples)
      output: temp("combined_19.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_19.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chr19 -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_X:
      input: expand("{x}_X.g.vcf",x=samples)
      output: temp("combined_X.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_X.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chrX -nt {threads}"+var
        print(cmd)
        shell(cmd)
        
  rule gatk_genotype_gvcfs_Y:
      input: expand("{x}_Y.g.vcf",x=samples)
      output: temp("combined_Y.vcf")
      params: gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:genGvcf"
      threads: 4
      run:
        fl=os.popen("ls *_Y.g.vcf").read().split()      
        var=" --variant "+" --variant ".join(fl)
        cmd="module load GATK/3.8-1; java -Xmx96g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T GenotypeGVCFs -R {params.genome} --disable_auto_index_creation_and_locking_when_reading_rods --annotation InbreedingCoeff --annotation FisherStrand --annotation QualByDepth --annotation ChromosomeCounts  --dbsnp {params.snpsites} -o {output} -L chrY -nt {threads}"+var
        print(cmd)
        shell(cmd)