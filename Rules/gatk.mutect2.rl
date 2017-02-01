if config['project']['annotation'] == "hg19":
  rule gatk_mutect2_1:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_1.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_1"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 1 -o {output.vcf}"

  rule gatk_mutect2_2:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_2.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_2"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 2 -o {output.vcf}"

  rule gatk_mutect2_3:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_3.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_3"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 3 -o {output.vcf}"

  rule gatk_mutect2_4:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_4.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_4"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 4 -o {output.vcf}"

  rule gatk_mutect2_5:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_5.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_5"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 5 -o {output.vcf}"

  rule gatk_mutect2_6:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_6.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_6"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 6 -o {output.vcf}"

  rule gatk_mutect2_7:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_7.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_7"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 7 -o {output.vcf}"

  rule gatk_mutect2_8:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_8.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_8"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 8 -o {output.vcf}"

  rule gatk_mutect2_9:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_9.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_9"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 9 -o {output.vcf}"

  rule gatk_mutect2_10:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_10.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_10"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 10 -o {output.vcf}"

  rule gatk_mutect2_11:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_11.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_11"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 11 -o {output.vcf}"

  rule gatk_mutect2_12:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_12.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_12"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 12 -o {output.vcf}"

  rule gatk_mutect2_13:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_13.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_13"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 13 -o {output.vcf}"

  rule gatk_mutect2_14:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_14.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_14"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 14 -o {output.vcf}"

  rule gatk_mutect2_15:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_15.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_15"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 15 -o {output.vcf}"

  rule gatk_mutect2_16:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_16.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_16"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 16 -o {output.vcf}"

  rule gatk_mutect2_17:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_17.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_17"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 17 -o {output.vcf}"

  rule gatk_mutect2_18:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_18.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_18"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 18 -o {output.vcf}"

  rule gatk_mutect2_19:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_19.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_19"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 19 -o {output.vcf}"

  rule gatk_mutect2_20:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_20.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_20"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 20 -o {output.vcf}"

  rule gatk_mutect2_21:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_21.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_21"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 21 -o {output.vcf}"

  rule gatk_mutect2_22:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_22.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_22"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L 22 -o {output.vcf}"

  rule gatk_mutect2_X:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_X.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_X"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L X -o {output.vcf}"

  rule gatk_mutect2_Y:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_Y.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_Y"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L Y -o {output.vcf}"

else:

  rule gatk_mutect2_1:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr1.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr1"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr1 -o {output.vcf}"

  rule gatk_mutect2_2:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr2.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr2"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr2 -o {output.vcf}"

  rule gatk_mutect2_3:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr3.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr3"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr3 -o {output.vcf}"

  rule gatk_mutect2_4:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr4.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr4"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr4 -o {output.vcf}"

  rule gatk_mutect2_5:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr5.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr5"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr5 -o {output.vcf}"

  rule gatk_mutect2_6:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr6.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr6"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr6 -o {output.vcf}"

  rule gatk_mutect2_7:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr7.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr7"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr7 -o {output.vcf}"

  rule gatk_mutect2_8:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr8.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr8"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr8 -o {output.vcf}"

  rule gatk_mutect2_9:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr9.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr9"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr9 -o {output.vcf}"

  rule gatk_mutect2_10:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr10.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr10"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr10 -o {output.vcf}"

  rule gatk_mutect2_11:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr11.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr11"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr11 -o {output.vcf}"

  rule gatk_mutect2_12:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr12.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr12"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr12 -o {output.vcf}"

  rule gatk_mutect2_13:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr13.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr13"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr13 -o {output.vcf}"

  rule gatk_mutect2_14:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr14.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr14"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr14 -o {output.vcf}"

  rule gatk_mutect2_15:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr15.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr15"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr15 -o {output.vcf}"

  rule gatk_mutect2_16:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr16.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr16"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr16 -o {output.vcf}"

  rule gatk_mutect2_17:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr17.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr17"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr17 -o {output.vcf}"

  rule gatk_mutect2_18:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr18.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr18"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr18 -o {output.vcf}"

  rule gatk_mutect2_19:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chr19.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chr19"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chr19 -o {output.vcf}"

  rule gatk_mutect2_X:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chrX.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chrX"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chrX -o {output.vcf}"

  rule gatk_mutect2_Y:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           targets="exome_targets.bed"
    output: vcf=config['project']['workpath']+"/mutect2_out/chrom_files/{x}_chrY.vcf"
    params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],knowns=config['references'][pfamily]['MUTECTVARIANTS'],targets="exome_targets.bed",rname="pl:mutect2_chrY"
    threads: 1
    shell: "module load GATK/3.7; GATK -m 120G MuTect2 -R {params.genome} -I:tumor {input.tumor} -I:normal {input.normal} --read_filter BadCigar --annotation Coverage -A FisherStrand -A HaplotypeScore -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest {params.knowns} -L chrY -o {output.vcf}"
## sed -i 's/NORMAL/{params.normalsample}/g' {output.vcfRename}; sed -i 's/TUMOR/{params.tumorsample}/g' {output.vcfRename}