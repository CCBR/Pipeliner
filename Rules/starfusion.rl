rule starfusion:
     input:  file1="{x}.R1.trimmed.fastq.gz",file2="{x}.R2.trimmed.fastq.gz"
     output: "starfusion/{x}/star-fusion.fusion_candidates.final"
     params: rname='starfusion',sample="{x}",starlib=config['references'][pfamily]['STARFUSIONLIB']
     shell: "module load STAR-Fusion/1.2.0; module load STAR/2.5.2a; mkdir -p starfusion/{params.sample}; STAR-Fusion --genome_lib_dir {params.starlib} --left_fq {input.file1} --right_fq {input.file2} --output_dir starfusion/{params.sample}"