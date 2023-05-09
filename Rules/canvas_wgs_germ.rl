rule canvas_wgs_germ:
    input: bam=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam",
           vcf="sample_vcfs/{x}.sample.vcf",
           index1="{x}.recal.bam.bai",
           index2="{x}.recal.bai"
    output: vcf="canvas_out/{x}/CNV.vcf.gz",
    params: genome=config['references'][pfamily]['CANVASGENOME'],ploidy=config['references'][pfamily]['CANVASPLOIDY'],sample=lambda wildcards: config['project']['units'][wildcards.x],kmer=config['references'][pfamily]['CANVASKMER'],filter=config['references'][pfamily]['CANVASFILTER'],balleles=config['references'][pfamily]['CANVASBALLELES'],rname="pl:canvas"
    shell: "mkdir -p canvas_out; mkdir -p canvas_out/{params.sample}; cp {params.ploidy} canvas_out/{params.sample}/; sed -i 's/SAMPLENAME/{params.sample}/g' canvas_out/{params.sample}/ploidy.vcf; export COMPlus_gcAllowVeryLargeObjects=1; module load Canvas/1.38; Canvas.dll Germline-WGS -b {input.bam} -n {params.sample} -o canvas_out/{params.sample} -r {params.kmer} --ploidy-vcf=canvas_out/{params.sample}/ploidy.vcf -g {params.genome} -f {params.filter} --population-b-allele-vcf={params.balleles}"