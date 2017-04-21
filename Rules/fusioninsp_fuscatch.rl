rule fusioninsp_fuscatch:
     input:  fusions="fusioncatcher/{x}/final-list_candidate-fusion-genes.GRCh37.txt",file1="{x}.R1.trimmed.fastq.gz",file2="{x}.R2.trimmed.fastq.gz"
     output: "fusioncatcher/fusioninspector/{x}/{x}.fusion_predictions.final"
     params: rname='fusioninsp',sample="{x}",starlib=config['references'][pfamily]['STARFUSIONLIB']
     threads: 8
     shell: module load fusioninspector/1.0.1; mkdir fusioncatcher/fusioninspector/{params.sample}; perl Rules/make_fusioninspector_input.pl {params.sample}; FusionInspector --fusions fusioncatcher/fusioninspector/{params.sample}/{params.sample}_fusionInspector.input --genome_lib {params.starlib} --left-fq {input.file1} --right-fq {input.file2} --out-dir fusioncatcher/fusioninspector/{params.sample} --out_prefix {params.sample} --prep_for_IGV