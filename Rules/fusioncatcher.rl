rule fusioncatcher:
     input:  file1="{x}.R1."+config['project']['filetype'],file2="{x}.R2."+config['project']['filetype']
     output: "fusioncatcher/{x}/final-list_candidate-fusion-genes.txt"
     params: data=config['references'][pfamily]['FUSCATCHDAT'],configfile=config['references'][pfamily]['FUSCATCHDAT'],rname='fusioncatcher',sample="{x}"
     threads: 16
     shell: "module load fusioncatcher/1.10; mkdir -p fusioncatcher/{params.sample}; fusioncatcher -i {input.file1},{input.file2} -o fusioncatcher/{params.sample} --threads {threads}"
