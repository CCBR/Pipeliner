if config['project']['annotation'] == "hg19":

  rule gatk_haplotype_caller_1:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_1.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 1"

  rule gatk_haplotype_caller_2:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_2.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 2"

  rule gatk_haplotype_caller_3:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_3.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 3"

  rule gatk_haplotype_caller_4:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_4.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 4"

  rule gatk_haplotype_caller_5:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_5.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 5"

  rule gatk_haplotype_caller_6:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_6.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 6"

  rule gatk_haplotype_caller_7:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_7.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 7"

  rule gatk_haplotype_caller_8:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_8.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 8"

  rule gatk_haplotype_caller_9:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_9.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 9"

  rule gatk_haplotype_caller_10:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_10.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 10"

  rule gatk_haplotype_caller_11:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_11.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 11"

  rule gatk_haplotype_caller_12:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_12.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 12"

  rule gatk_haplotype_caller_13:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_13.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 13"

  rule gatk_haplotype_caller_14:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_14.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 14"

  rule gatk_haplotype_caller_15:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_15.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 15"

  rule gatk_haplotype_caller_16:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_16.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 16"

  rule gatk_haplotype_caller_17:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_17.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 17"

  rule gatk_haplotype_caller_18:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_18.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 18"

  rule gatk_haplotype_caller_19:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_19.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 19"

  rule gatk_haplotype_caller_20:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_20.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 20"

  rule gatk_haplotype_caller_21:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_21.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 21"

  rule gatk_haplotype_caller_22:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_22.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L 22"

  rule gatk_haplotype_caller_X:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_X.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L X"

  rule gatk_haplotype_caller_Y:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_Y.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L Y"

  rule gatk_haplotype_caller_MT:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_MT.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -ploidy 1 -L MT"

elif config['project']['annotation'] == "hg38":

  rule gatk_haplotype_caller_1:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_1.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr1"

  rule gatk_haplotype_caller_2:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_2.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr2"

  rule gatk_haplotype_caller_3:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_3.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr3"

  rule gatk_haplotype_caller_4:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_4.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr4"

  rule gatk_haplotype_caller_5:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_5.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr5"

  rule gatk_haplotype_caller_6:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_6.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr6"

  rule gatk_haplotype_caller_7:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_7.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr7"

  rule gatk_haplotype_caller_8:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_8.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr8"

  rule gatk_haplotype_caller_9:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_9.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr9"

  rule gatk_haplotype_caller_10:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_10.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr10"

  rule gatk_haplotype_caller_11:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_11.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr11"

  rule gatk_haplotype_caller_12:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_12.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr12"

  rule gatk_haplotype_caller_13:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_13.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr13"

  rule gatk_haplotype_caller_14:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_14.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr14"

  rule gatk_haplotype_caller_15:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_15.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr15"

  rule gatk_haplotype_caller_16:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_16.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr16"

  rule gatk_haplotype_caller_17:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_17.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr17"

  rule gatk_haplotype_caller_18:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_18.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr18"

  rule gatk_haplotype_caller_19:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_19.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr19"

  rule gatk_haplotype_caller_20:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_20.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr20"

  rule gatk_haplotype_caller_21:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_21.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr21"

  rule gatk_haplotype_caller_22:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_22.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr22"

  rule gatk_haplotype_caller_X:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_X.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chrX"

  rule gatk_haplotype_caller_Y:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_Y.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chrY"

  rule gatk_haplotype_caller_MT:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_MT.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -ploidy 1 -L chrM"

elif config['project']['annotation'] == "mm10":

  rule gatk_haplotype_caller_1:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_1.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr1"

  rule gatk_haplotype_caller_2:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_2.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr2"

  rule gatk_haplotype_caller_3:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_3.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr3"

  rule gatk_haplotype_caller_4:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_4.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr4"

  rule gatk_haplotype_caller_5:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_5.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr5"

  rule gatk_haplotype_caller_6:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_6.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr6"

  rule gatk_haplotype_caller_7:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_7.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr7"

  rule gatk_haplotype_caller_8:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_8.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr8"

  rule gatk_haplotype_caller_9:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_9.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr9"

  rule gatk_haplotype_caller_10:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_10.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr10"

  rule gatk_haplotype_caller_11:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_11.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr11"

  rule gatk_haplotype_caller_12:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_12.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr12"

  rule gatk_haplotype_caller_13:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_13.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr13"

  rule gatk_haplotype_caller_14:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_14.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr14"

  rule gatk_haplotype_caller_15:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_15.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr15"

  rule gatk_haplotype_caller_16:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_16.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr16"

  rule gatk_haplotype_caller_17:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_17.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr17"

  rule gatk_haplotype_caller_18:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_18.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr18"

  rule gatk_haplotype_caller_19:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_19.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chr19"

  rule gatk_haplotype_caller_X:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_X.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chrX"

  rule gatk_haplotype_caller_Y:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_Y.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -L chrY"

  rule gatk_haplotype_caller_MT:
     input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam", index=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bai"
     output:temp("{x}_MT.g.vcf")
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:hapcall"
     threads: 4
     shell: "module load GATK/3.8-1; java -Xmx64g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T HaplotypeCaller -R {params.genome} -I {input.bam} --read_filter BadCigar --use_jdk_inflater --use_jdk_deflater --emitRefConfidence GVCF --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest --variant_index_type LINEAR --dbsnp {params.snpsites} --variant_index_parameter 128000 -nct {threads} -o {output[0]} -ploidy 1 -L chrM"
