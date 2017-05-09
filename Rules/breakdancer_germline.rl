rule breakdancer_germline:
     input: expand(config['project']['workpath']+"/{x}.recal.bam", x=samples)
     output: "breakdancer_out/file.ctx"
     params: genome=config['references'][pfamily]['GENOME'],rname="breakdancer"
     threads: 8
     shell: "mkdir -p breakdancer_out; module load breakdancer/1.4.5; cd breakdancer_out; bam2cfg.pl -g -h {input} > config_file.cfg; breakdancer-max config_file.cfg -g support.bed > file.ctx"