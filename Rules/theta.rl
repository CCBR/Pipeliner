rule theta:
       input:  normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
               tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
               calls="cnvkit_out/{x}_calls.cns",
               targets="exome_targets.bed"
       output: infile="{x}.input",
               normalcon=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".config",
               tumorcon=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".config",
               normalsnpfile=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".withCounts",
               tumorsnpfile=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".withCounts",
               graph=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".input.n2.graph.pdf"
       params: normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],genome=config['references'][pfamily]['GENOME'],exons=config['references'][pfamily]['EXONS'],sites=config['references'][pfamily]['THETASNPS'],rname="pl:theta"
       threads: 4
       shell:  "module load samtools; module load theta; CreateExomeInput -s cnvkit_out/{params.normalsamples}+{params.tumorsample}_cnvkit/{params.tumorsample}.cns -t {input.tumor} -n {input.normal} --OUTPUT_PREFIX {params.tumorsample} --FA {params.genome} --EXON_FILE {params.exons} --DIR theta_out; perl Scripts/make_theta_config.pl {params.normalsample} {params.tumorsample} {params.sites}; java -Xmx24g -jar $THETA_JARPATH/getAlleleCounts.jar theta_out/{output.tumorcon}; java -Xmx24g -jar $THETA_JARPATH/getAlleleCounts.jar theta_out/{output.normalcon}; RunTHetA theta_out/{output.infile} --TUMOR_SNP theta_out/{output.tumorsnpfile} --NORMAL_FILE theta_out/{output.normalsnpfile} --NUM_PROCESSES {threads} --OUTPUT_PREFIX {params.tumorsample} --DIR theta_out --BAF"