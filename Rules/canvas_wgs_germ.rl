rule canvas_wgs_germ:
    input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam",
           vcf="sample_vcfs/{x}.sample.vcf",
    output: vcf="canvas_out/{x}/CNV.vcf.gz",
    params: genome=config['references'][pfamily]['CANVASGENOME'],sample=lambda wildcards: config['project']['units'][wildcards.x],kmer=config['references'][pfamily]['CANVASKMER'],filter=config['references'][pfamily]['CANVASFILTER'],rname="pl:canvas"
    shell: "mkdir -p canvas_out; mkdir -p canvas_out/{params.sample}; cp {params.sample}.recal.bai {params.sample}.recal.bam.bai; export COMPlus_gcAllowVeryLargeObjects=1; module load Canvas/1.25; Canvas.dll Germline-WGS -b {input.bam} -o canvas_out/{params.sample} -r {params.kmer} -g {params.genome} -f {params.filter} --sample-name={params.sample} --sample-b-allele-vcf={input.vcf}"