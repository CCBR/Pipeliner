rule oncofuse_fuscatch:
     input:  fusions="fusioncatcher/fusioninspector/{x}/{x}.fusion_predictions.final"
     output: infile="fusioncatcher/oncofuse/{x}/{x}.oncofuse.input",outfile="fusioncatcher/oncofuse/{x}/{x}.oncofuse.output"
     params: rname='oncofuse',sample="{x}",oncofuse=config['bin'][pfamily]['ONCOFUSE'],build=config['references'][pfamily]['ONCOFUSEBUILD']
     shell: "mkdir -p fusioncatcher/oncofuse/{params.sample}; module load perl/5.18.4; perl Scripts/make_oncofuse_input.pl {params.sample} fusioncatcher; {params.oncofuse} -p 1 -a {params.build} {output.infile} coord - {output.outfile}"