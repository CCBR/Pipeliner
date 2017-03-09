rule purity:
    input: normal=lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam",
           tumor=lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam",
           calls="cnvkit_out/{x}_calls.cns"
    output: gainloss="cnvkit_out/{x}_gainloss.tsv"
    params: tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],genome=config['references'][pfamily]['GENOME'],exons=config['references'][pfamily]['EXONS'],snps=config['references'][pfamily]['THETASNPS'],rname="pl:purity"
    shell: "module load samtools; module load theta; CreateExomeInput -s cnvkit_out/{params.tumorsample}.cns -t {input.tumor} -n {input.normal} --OUTPUT_PREFIX {params.tumorsample} --FA {params.genome} --EXON_FILE {params.exons}"