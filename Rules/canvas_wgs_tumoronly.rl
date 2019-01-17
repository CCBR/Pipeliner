rule canvas_wgs_tumoronly:
    input: tumor=lambda wildcards: config['project']['units'][wildcards.x]+".recal.bam",
           somaticvcf=config['project']['workpath']+"/mutect2_out/{x}.FINALmutect2.vcf",
           index1="{x}.recal.bam.bai",
           index2="{x}.recal.bai"
    output: vcf="canvas_out/{x}/CNV.vcf.gz",
    params: genome=config['references'][pfamily]['CANVASGENOME'],tumorsample=lambda wildcards: config['project']['units'][wildcards.x],kmer=config['references'][pfamily]['CANVASKMER'],filter=config['references'][pfamily]['CANVASFILTER'],balleles=config['references'][pfamily]['CANVASBALLELES'],ploidyvcf=config['references'][pfamily]['CANVASPLOIDY'],rname="pl:canvas"
    shell: "mkdir -p canvas_out; mkdir -p canvas_out/{params.tumorsample}; cp {params.ploidyvcf} canvas_out/{params.tumorsample}/; sed -i 's/SAMPLENAME/{params.tumorsample}/g' canvas_out/{params.tumorsample}/ploidy.vcf; export COMPlus_gcAllowVeryLargeObjects=1; module load Canvas/1.38; Canvas.dll Somatic-WGS -b {input.tumor} -n {params.tumorsample} -o canvas_out/{params.tumorsample} -r {params.kmer} -g {params.genome} -f {params.filter} --population-b-allele-vcf={params.balleles} --somatic-vcf={input.somaticvcf} --ploidy-vcf=canvas_out/{params.tumorsample}/ploidy.vcf"