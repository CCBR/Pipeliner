rule fusioninsp:
     input:  starfusion="starfusion/{x}/star-fusion.fusion_predictions.tsv",fusioncatcher="fusioncatcher/fusioninspector/{x}/{x}_fusionInspector.input",file1="{x}.R1.trimmed.fastq.gz",file2="{x}.R2.trimmed.fastq.gz"
     output: "fusioninspector/{x}/{x}.fusion_predictions.final"
     params: rname='fusioninsp',sample="{x}",ref=config['project']['annotation'],starlib=config['references'][pfamily]['STARFUSIONLIB']
     shell: "module load fusioninspector/1.1.0; module load STAR/2.7.0f; mkdir -p fusioninspector/{params.sample}; module load perl; FusionInspector --fusions {input.starfusion},{input.fusioncatcher} --genome_lib {params.starlib} --left_fq {input.file1} --right_fq {input.file2} --out_dir fusioninspector/{params.sample} --out_prefix {params.sample} --prep_for_IGV --CPU 16 --cleanup --annotate --examine_coding_effect"
