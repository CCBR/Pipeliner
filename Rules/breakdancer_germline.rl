rule breakdancer_germline:
     input: expand("{x}.recal.bam", x=samples)
     output: "breakdancer_out/file.ctx"
     params: genome=config['references'][pfamily]['GENOME'],rname="breakdancer"
     threads: 8
     shell: "mkdir -p breakdancer_out; module load breakdancer/1.4.5; bam2cfg.pl -g -h {input} > breakdancer_out/config_file.cfg; cd breakdancer_out; breakdancer-max config_file.cfg -g support.bed > file.ctx"