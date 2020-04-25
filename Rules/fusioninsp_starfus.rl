rule fusioninsp_starfus:
     input:  fusions="starfusion/{x}/star-fusion.fusion_predictions.tsv",file1="{x}.R1.trimmed.fastq.gz",file2="{x}.R2.trimmed.fastq.gz"
     output: "starfusion/fusioninspector/{x}/{x}.fusion_predictions.final"
     params: rname='fusioninsp',sample="{x}",starlib=config['references'][pfamily]['STARFUSIONLIB']
     shell: """
numfusions=$(grep -v ^# {input.fusions} | wc -l)
if (( $numfusions > 0 )); then
     module load fusioninspector/1.1.0; module load STAR/2.7.0f; mkdir -p starfusion/fusioninspector/{params.sample}; module load perl; FusionInspector --fusions {input.fusions} --genome_lib {params.starlib} --left_fq {input.file1} --right_fq {input.file2} --out_dir starfusion/fusioninspector/{params.sample} --out_prefix {params.sample} --prep_for_IGV
else
     touch {output};
fi
"""
