rule trimmomatic:
    input:  "{x}.R1."+config['project']['filetype'],
            "{x}.R2."+config['project']['filetype']
    output: "{x}.R1.trimmed.fastq.gz",
            "{x}.R1.trimmed.unpair.fastq.gz",
            "{x}.R2.trimmed.fastq.gz",
            "{x}.R2.trimmed.unpair.fastq.gz"
    params: trimmomatic=config['bin']['TRIMMOMATIC'],
            adapterfile=config['references']['trimmomatic.adapters'],rname="pl:trimmomatic"
    threads: 4
    shell:  """
            {params.trimmomatic} PE -threads {threads} -phred33 {input[0]} {input[1]} {output} ILLUMINACLIP:{params.adapterfile}:3:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:20
           """

# MAXINFO:50:0.8     



