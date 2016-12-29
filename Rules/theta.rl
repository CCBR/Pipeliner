rule theta:
       input:  normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
               tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
               calls="cnvkit_out/{x}_calls.cns",
               targets="exome_targets.bed",
               dir="theta_out"
       output: infile="theta_out/{x}.input"
       params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],genome=config['references'][pfamily]['GENOME'],exons=config['references'][pfamily]['EXONS'],sites=config['references'][pfamily]['THETASNPS'],rname="pl:theta"
       threads: 8
       shell:  "module load samtools; module load theta; CreateExomeInput -s cnvkit_out/{params.normalsample}+{params.tumorsample}_cnvkit/{params.tumorsample}.cns -t {input.tumor} -n {input.normal} --OUTPUT_PREFIX {params.normalsample}+{params.tumorsample} --FA {params.genome} --EXON_FILE {params.exons} --DIR=theta_out/; perl Scripts/make_theta_config.pl {params.normalsample} {params.tumorsample} {params.sites}; java -Xmx24g -jar $THETA_JARPATH/getAlleleCounts.jar theta_out/{params.tumorsample}.config; java -Xmx24g -jar $THETA_JARPATH/getAlleleCounts.jar theta_out/{params.normalsample}.config; RunTHetA theta_out/{params.normalsample}+{params.tumorsample}.input --TUMOR_SNP theta_out/{params.tumorsample}.withCounts --NORMAL_FILE theta_out/{params.normalsample}.withCounts --NUM_PROCESSES {threads} --OUTPUT_PREFIX {params.tumorsample} --DIR theta_out/ --BAF"