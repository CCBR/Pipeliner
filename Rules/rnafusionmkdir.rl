rule rnafusionmkdir:
        input:  expand("{x}.R1.fastq.gz", x=samples),
        output: fusioncatcher=config['project']['workpath']+"/fusioncatcher",
                fusioninspcatch=config['project']['workpath']+"/fusioncatcher/fusioninspector",
                oncofusecatch=config['project']['workpath']+"/fusioncatcher/oncofuse",
                starfusion=config['project']['workpath']+"/starfusion",
                fusioninspstar=config['project']['workpath']+"/starfusion/fusioninspector",
                oncofusestar=config['project']['workpath']+"/starfusion/oncofuse",
                expression=config['project']['workpath']+"/expression",
        params: rname="pl:mkdir"
        shell:  "mkdir -p {output.fusioncatcher}; mkdir -p {output.starfusion}; mkdir -p {output.fusioninspcatch}; mkdir -p {output.oncofusecatch}; mkdir -p {output.fusioninspstar}; mkdir -p {output.oncofusestar}; mkdir -p {output.expression}"