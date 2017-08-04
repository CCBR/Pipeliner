rule trimmomatic:
    input:  "{x}.R1."+config['project']['filetype'],
            "{x}.R2."+config['project']['filetype']
    output: one=temp("{x}.R1.trimmed.fastq.gz"),
            two=temp("{x}.R1.trimmed.unpair.fastq.gz"),
            three=temp("{x}.R2.trimmed.fastq.gz"),
            four=temp("{x}.R2.trimmed.unpair.fastq.gz"),
            err="QC/{x}_run_trimmomatic.err"
    params: trimmomatic=config['bin'][pfamily]['TRIMMOMATIC'],
            adapterfile=config['references'][pfamily]['trimmomatic.adapters'],rname="pl:trimmomatic"
    threads: 32
    shell:  """
            {params.trimmomatic} PE -threads {threads} -phred33 {input[0]} {input[1]} {output.one} {output.two} {output.three} {output.four} ILLUMINACLIP:{params.adapterfile}:3:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:20 2> {output.err}
           """

# MAXINFO:50:0.8     



