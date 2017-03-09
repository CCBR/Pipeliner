rule novocraft_novoalign:
     input:  "{x}.R1."+config['project']['filetype'],
             "{x}.R2."+config['project']['filetype']
     output: bam = temp("{x}.bam"),
             qcal = "{x}.qcal"
     params: partition="norm", mem="16g", time="24:00:00",
             adapter1=config['references'][pfamily]['ADAPTER1'],
             adapter2=config['references'][pfamily]['ADAPTER2'],
             novoindex=config['references'][pfamily]['NOVOINDEX'],
             qcal=config['bin'][pfamily]['QCAL'],
             sam=config['bin'][pfamily]['SAMTOOLS'],rname="pl:novoalign"
#     message: "Executing NovoAlign Rule."
     threads: 32
     version: "1.0"
     shell: "module load novocraft/3.02.10; novoalign -c {threads} -d {params.novoindex} -a {params.adapter1} {params.adapter2} -H 10 -k -K {output.qcal} -i PE 180,50 -o SAM -f {input} | {params.sam} view -bS - > {output.bam};"
#     shell: "module load novocraft/3.02.10; novoalign -c {threads} -d {params.novoindex} -a {params.adapter1} {params.adapter2} -H 10 -k -K {output.qcal} -i PE 180,50 -o SAM -f {input} > {output.sam};"

#{params.qcal} {output.qcal} {output.qcal}.qcalreport.pdf;
