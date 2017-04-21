rule rnafusionmkdir:
        input:  "{x}.R1."+config['project']['filetype'],
                "{x}.R2."+config['project']['filetype']
        output: fusioncatcher=config['project']['workpath']+"/fusioncatcher",
                fusioninspcatch=config['project']['workpath']+"/fusioncatcher/fusioninspector",
                oncofusecatch=config['project']['workpath']+"/fusioncatcher/oncofuse",
                starfusion=config['project']['workpath']+"/starfusion",
                fusioninspstar=config['project']['workpath']+"/starfusion/fusioninspector",
                oncofusestar=config['project']['workpath']+"/starfusion/oncofuse",
        params: rname="pl:mkdir"
        shell:  "mkdir -p {output.fusioncatcher}; mkdir -p {output.starfusion}; mkdir -p {output.fusioninspcatch}; mkdir -p {output.oncofusecatch}; mkdir -p {fusioninspstar}; mkdir -p {oncofusestar}