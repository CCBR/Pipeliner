if config['project']['annotation'] == "hg19":
 rule freec_wgs_somatic_pass2:
     input: newnormbam=temp("freec_out/"+lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
            newtumorbam=temp("freec_out/"+lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
            index1=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bai",
            index2=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bai",
            fit="sequenza_out/{x}/{x}"+"_alternative_solutions.txt",
     output: cnvs="freec_out/pass2/{x}/"+lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam_CNVs",
             newnormbam=temp("freec_out/"+lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam"),
             newtumorbam=temp("freec_out/"+lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam"),
             normalpileup="freec_out/pass2/{x}/"+lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam_minipileup.pileup",
             tumorpileup="freec_out/pass2/{x}/"+lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam_minipileup.pileup",
     params: dir=config['project']['workpath'],fasta=config['references'][pfamily]['FREECFASTA'],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],lengths=config['references'][pfamily]['FREECLENGTHS'],chroms=config['references'][pfamily]['FREECCHROMS'],pile=config['references'][pfamily]['FREECPILEUP'],snps=config['references'][pfamily]['FREECSNPS'],rname="pl:freec"
     shell: "module load samtools/1.9
             mkdir -p freec_out
             mkdir -p freec_out/pass2
             mkdir -p freec_out/pass2/{params.normalsample}+{params.tumorsample}
             perl Scripts/make_freec_pass2_wgs_tn_config.pl {params.dir}/freec_out/pass2/{params.normalsample}+{params.tumorsample} {params.lengths} {params.chroms} {params.dir}/freec_out/{input.newtumorbam} {params.dir}/freec_out/{input.newnormalbam} {params.pile} {params.fasta} {params.snps} {input.fit}
             module load freec/11.5
             module load bedtools/2.27.1
             freec -conf {params.dir}/freec_out/pass2/{params.normalsample}+{params.tumorsample}/freec_wgs_config.txt"

else:
 rule freec_wgs_somatic_pass2:
     input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
            tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
            index1=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bai",
            index2=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bai",
            fit="sequenza_out/{x}/{x}"+"_alternative_solutions.txt",
     output: cnvs="freec_out/pass2/{x}/"+lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam_CNVs",
             normalpileup="freec_out/pass2/{x}/"+lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam_minipileup.pileup",
             tumorpileup="freec_out/pass2/{x}/"+lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam_minipileup.pileup",
     params: dir=config['project']['workpath'],fasta=config['references'][pfamily]['FREECFASTA'],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],lengths=config['references'][pfamily]['FREECLENGTHS'],chroms=config['references'][pfamily]['FREECCHROMS'],pile=config['references'][pfamily]['FREECPILEUP'],snps=config['references'][pfamily]['FREECSNPS'],rname="pl:freec"
     shell: "mkdir -p freec_out
             mkdir -p freec_out/pass2
             mkdir -p freec_out/pass2/{params.normalsample}+{params.tumorsample}
             perl Scripts/make_freec_pass2_wgs_tn_config.pl {params.dir}/freec_out/pass2/{params.normalsample}+{params.tumorsample} {params.lengths} {params.chroms} {params.dir}/{input.tumor} {params.dir}/{input.normal} {params.pile} {params.fasta} {params.snps} {input.fit}
             module load freec/11.5
             module load samtools/1.9
             module load bedtools/2.27.1
             freec -conf {params.dir}/freec_out/pass2/{params.normalsample}+{params.tumorsample}/freec_wgs_config.txt"