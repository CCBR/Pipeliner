rule rule oncofuse_starfus:
     input:  fusions="starfusion/fusioninspector/{x}/{x}.fusion_predictions.final"
     output: infile="starfusion/oncofuse/{x}/{x}.oncofuse.input",outfile="starfusion/oncofuse/{x}/{x}.oncofuse.output"
     params: rname='oncofuse',sample="{x}",oncofuse=config['bin'][pfamily]['ONCOFUSE'],build=config['references'][pfamily]['ONCOFUSEBUILD']
     shell: "mkdir -p starfusion/oncofuse/{params.sample}; module load perl/5.18.4; perl Scripts/make_oncofuse_input.pl {params.sample} starfusion; {params.oncofuse} -p 1 -a {params.build} {output.infile} coord - {output.outfile}"