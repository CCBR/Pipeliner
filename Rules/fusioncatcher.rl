rule fusioncatcher:
     input:  file1="{x}.R1."+config['project']['filetype'],file2="{x}.R2."+config['project']['filetype']
     output: "fusioncatcher/{x}/final-list_candidate-fusion-genes.GRCh37.txt"
     params: rname='fusioncatcher',sample="{x}"
     threads: 16
     shell: "module load fusioncatcher/0.99.7b; mkdir fusioncatcher/{params.sample}; fusioncatcher -i {input.file1},{input.file2} -o fusioncatcher/{params.sample} --threads {threads}"