if config['project']['annotation'] == "hg19":
  rule gatk_mutect2_1:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_1.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_1"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 1 -o {output.vcf}"

  rule gatk_mutect2_2:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_2.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_2"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 2 -o {output.vcf}"

  rule gatk_mutect2_3:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_3.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_3"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 3 -o {output.vcf}"

  rule gatk_mutect2_4:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_4.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_4"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 4 -o {output.vcf}"

  rule gatk_mutect2_5:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_5.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_5"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 5 -o {output.vcf}"

  rule gatk_mutect2_6:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_6.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_6"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 6 -o {output.vcf}"

  rule gatk_mutect2_7:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_7.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_7"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 7 -o {output.vcf}"

  rule gatk_mutect2_8:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_8.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_8"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 8 -o {output.vcf}"

  rule gatk_mutect2_9:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_9.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_9"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 9 -o {output.vcf}"

  rule gatk_mutect2_10:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_10.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_10"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 10 -o {output.vcf}"

  rule gatk_mutect2_11:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_11.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_11"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 11 -o {output.vcf}"

  rule gatk_mutect2_12:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_12.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_12"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 12 -o {output.vcf}"

  rule gatk_mutect2_13:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_13.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_13"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 13 -o {output.vcf}"

  rule gatk_mutect2_14:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_14.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_14"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 14 -o {output.vcf}"

  rule gatk_mutect2_15:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_15.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_15"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 15 -o {output.vcf}"

  rule gatk_mutect2_16:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_16.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_16"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 16 -o {output.vcf}"

  rule gatk_mutect2_17:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_17.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_17"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 17 -o {output.vcf}"

  rule gatk_mutect2_18:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_18.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_18"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 18 -o {output.vcf}"

  rule gatk_mutect2_19:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_19.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_19"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 19 -o {output.vcf}"

  rule gatk_mutect2_20:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_20.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_20"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 20 -o {output.vcf}"

  rule gatk_mutect2_21:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_21.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_21"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 21 -o {output.vcf}"

  rule gatk_mutect2_22:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_22.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_22"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 22 -o {output.vcf}"

  rule gatk_mutect2_X:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_X.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_X"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L X -o {output.vcf}"

  rule gatk_mutect2_Y:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_Y.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_Y"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L Y -o {output.vcf}"

elif config['project']['annotation'] == "hg38":

  rule gatk_mutect2_1:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_1.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_1"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr1 -o {output.vcf}"

  rule gatk_mutect2_2:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_2.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_2"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr2 -o {output.vcf}"

  rule gatk_mutect2_3:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_3.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_3"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr3 -o {output.vcf}"

  rule gatk_mutect2_4:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_4.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_4"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr4 -o {output.vcf}"

  rule gatk_mutect2_5:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_5.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_5"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr5 -o {output.vcf}"

  rule gatk_mutect2_6:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_6.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_6"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr6 -o {output.vcf}"

  rule gatk_mutect2_7:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_7.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_7"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr7 -o {output.vcf}"

  rule gatk_mutect2_8:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_8.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_8"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr8 -o {output.vcf}"

  rule gatk_mutect2_9:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_9.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_9"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr9 -o {output.vcf}"

  rule gatk_mutect2_10:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_10.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_10"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr10 -o {output.vcf}"

  rule gatk_mutect2_11:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_11.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_11"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr11 -o {output.vcf}"

  rule gatk_mutect2_12:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_12.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_12"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr12 -o {output.vcf}"

  rule gatk_mutect2_13:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_13.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_13"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr13 -o {output.vcf}"

  rule gatk_mutect2_14:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_14.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_14"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr14 -o {output.vcf}"

  rule gatk_mutect2_15:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_15.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_15"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr15 -o {output.vcf}"

  rule gatk_mutect2_16:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_16.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_16"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr16 -o {output.vcf}"

  rule gatk_mutect2_17:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_17.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_17"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr17 -o {output.vcf}"

  rule gatk_mutect2_18:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_18.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_18"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr18 -o {output.vcf}"

  rule gatk_mutect2_19:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_19.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_19"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr19 -o {output.vcf}"

  rule gatk_mutect2_20:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_20.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_20"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr20 -o {output.vcf}"

  rule gatk_mutect2_21:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_21.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_21"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr21 -o {output.vcf}"

  rule gatk_mutect2_22:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_22.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_22"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr22 -o {output.vcf}"

  rule gatk_mutect2_X:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_X.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_X"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chrX -o {output.vcf}"

  rule gatk_mutect2_Y:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_Y.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_Y"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chrY -o {output.vcf}"

elif config['project']['annotation'] == "mm10":

  rule gatk_mutect2_1:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_1.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr1"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr1 -o {output.vcf}"

  rule gatk_mutect2_2:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_2.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr2"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr2 -o {output.vcf}"

  rule gatk_mutect2_3:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_3.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr3"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr3 -o {output.vcf}"

  rule gatk_mutect2_4:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_4.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr4"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr4 -o {output.vcf}"

  rule gatk_mutect2_5:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_5.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr5"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr5 -o {output.vcf}"

  rule gatk_mutect2_6:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_6.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr6"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr6 -o {output.vcf}"

  rule gatk_mutect2_7:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_7.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr7"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr7 -o {output.vcf}"

  rule gatk_mutect2_8:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_8.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr8"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr8 -o {output.vcf}"

  rule gatk_mutect2_9:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_9.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr9"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr9 -o {output.vcf}"

  rule gatk_mutect2_10:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_10.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr10"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr10 -o {output.vcf}"

  rule gatk_mutect2_11:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_11.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr11"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr11 -o {output.vcf}"

  rule gatk_mutect2_12:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_12.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr12"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr12 -o {output.vcf}"

  rule gatk_mutect2_13:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_13.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr13"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr13 -o {output.vcf}"

  rule gatk_mutect2_14:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_14.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr14"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr14 -o {output.vcf}"

  rule gatk_mutect2_15:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_15.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr15"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr15 -o {output.vcf}"

  rule gatk_mutect2_16:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_16.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr16"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr16 -o {output.vcf}"

  rule gatk_mutect2_17:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_17.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr17"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr17 -o {output.vcf}"

  rule gatk_mutect2_18:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_18.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr18"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr18 -o {output.vcf}"

  rule gatk_mutect2_19:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_19.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chr19"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr19 -o {output.vcf}"

  rule gatk_mutect2_X:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_X.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chrX"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chrX -o {output.vcf}"

  rule gatk_mutect2_Y:
    input: normal=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
           tumor=ancient(lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_Y.vcf"
    params: pon=config['references'][pfamily]['PON'],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets=ancient("exome_targets.bed"),rname="pl:mutect2_chrY"
    threads: 2
    shell: "module load GATK/3.7; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} -PON {params.pon} -G Standard --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chrY -o {output.vcf}"