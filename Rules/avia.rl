rule avia:
    input: config['project']['workpath']+"/variants.bed"
    output:config['project']['workpath']+"/full_annot.txt.zip"
    params: batch ="-l nodes=1:gpfs -q ccr",email=config['project']['username'],species=config['references'][pfamily]['AVIASET'],rname="avia"
    shell: """
         perl Scripts/avia.pl {input} {params.email} {params.species}

           """

