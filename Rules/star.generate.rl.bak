rule star_generate:
     input: 
     output: "dict.built"
     params: star=config['bin']['STAR'],
             genomeDir=config['references']['genomeDir'],
             gtf=config['references']['gencodeGtf'],
             fasta=config['references']['genomeRnaseq'],rname="pl:stargen"
     shell:  "module load STAR;{params.star} --runMode genomeGenerate --genomeDir {params.genomeDir} --sjdbOverhang 99 --genomeFastaFiles {params.genomeDir}/{params.fasta} --runThreadN 32 --sjdbGTFfile {params.genomeDir}/{params.gtf};touch dict.built"


