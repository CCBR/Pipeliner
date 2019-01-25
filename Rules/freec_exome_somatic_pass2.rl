if config['project']['annotation'] == "hg19":
 rule freec_exome_somatic_pass2:
     input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
            tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
            fit="sequenza_out/{x}"+"_alternative_solutions.txt",
            freectargets=config['project']['workpath']+"/freec_targets.bed",
     output: cnvs="freec_out/pass2/{x}.recal.bam_CNVs",
     params: dir=config['project']['workpath'],fasta=config['references'][pfamily]['FREECFASTA'],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],lengths=config['references'][pfamily]['FREECLENGTHS'],chroms=config['references'][pfamily]['FREECCHROMS'],pile=config['references'][pfamily]['FREECPILEUP'],snps=config['references'][pfamily]['FREECSNPS'],rname="pl:freec"
     shell: "mkdir -p freec_out/pass2; mkdir -p freec_out/pass2/{params.tumorsample}; perl Scripts/make_freec_pass2_exome_tn_config.pl {params.dir}/freec_out/pass2/{params.tumorsample} {params.lengths} {params.chroms} {params.dir}/freec_out/{input.tumor} {params.dir}/freec_out/{input.normal} {params.pile} {params.fasta} {params.snps} {input.freectargets} {input.fit}; module load samtools/1.9; module load freec/11.5; module load bedtools/2.27.1; freec -conf {params.dir}/freec_out/pass2/{params.tumorsample}/freec_exome_config.txt; mv {params.dir}/freec_out/pass2/{params.tumorsample}/{params.tumorsample}.recal.bam_CNVs {output.cnvs}"

else:
 rule freec_exome_somatic_pass2:
     input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
            tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
            index1=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bai",
            index2=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bai",
            fit="sequenza_out/{x}"+"_alternative_solutions.txt",
            freectargets=config['project']['workpath']+"/freec_targets.bed",
     output: cnvs="freec_out/pass2/{x}.recal.bam_CNVs",
     params: dir=config['project']['workpath'],fasta=config['references'][pfamily]['FREECFASTA'],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],lengths=config['references'][pfamily]['FREECLENGTHS'],chroms=config['references'][pfamily]['FREECCHROMS'],pile=config['references'][pfamily]['FREECPILEUP'],snps=config['references'][pfamily]['FREECSNPS'],rname="pl:freec"
     shell: "mkdir -p freec_out; mkdir -p freec_out/pass2; mkdir -p freec_out/pass2/{params.tumorsample}; perl Scripts/make_freec_pass2_exome_tn_config.pl {params.dir}/freec_out/pass2/{params.tumorsample} {params.lengths} {params.chroms} {params.dir}/{input.tumor} {params.dir}/{input.normal} {params.pile} {params.fasta} {params.snps} {input.freectargets} {input.fit}; module load freec/11.5; module load samtools/1.9; module load bedtools/2.27.1; freec -conf {params.dir}/freec_out/pass2/{params.tumorsample}/freec_exome_config.txt; mv {params.dir}/freec_out/pass2/{params.tumorsample}/{params.tumorsample}.recal.bam_CNVs {output.cnvs}"