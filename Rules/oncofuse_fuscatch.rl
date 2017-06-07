rule oncofuse_fuscatch:
     input:  fusions="fusioncatcher/fusioninspector/{x}/{x}.fusion_predictions.final"
     output: infile="fusioncatcher/oncofuse/{x}/{x}.oncofuse.input",outfile="fusioncatcher/oncofuse/{x}/{x}.oncofuse.output"
     params: rname='oncofuse',sample="{x}",oncofuse=config['bin'][pfamily]['ONCOFUSE']
     shell: "mkdir -p fusioncatcher/oncofuse/{params.sample}; perl Scripts/make_oncofuse_input.pl {params.sample} fusioncatcher; {params.oncofuse} {output.infile} coord - {output.outfile}"