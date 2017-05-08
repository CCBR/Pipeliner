rule manta_germline:
     input: expand("{x}.recal.bam", x=samples)
     output: "breakdancer_out/file.ctx"
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:manta"
     threads: 8
     run:
      fl=os.popen("ls *.recal.bam").read().split()      
      var=" --bam "+" --bam ".join(fl)
      cmd="mkdir -p breakdancer_out; module load breakdancer/1.4.5; bam2cfg.pl "+var"; breakdancer-max config_file.cfg > {output}"
      print(cmd)