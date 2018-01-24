rule manta_wgs:
     input: expand("{x}.recal.bam", x=samples)
     output: vcf="manta_out/results/variants/diploidSV.vcf.gz",
             dir="manta_out"
     params: gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:manta"
     threads: 8
     run:
      fl=os.popen("ls *.recal.bam").read().split()      
      var=" --bam "+" --bam ".join(fl)
      cmd="mkdir -p manta_out; module load python/2.7; module load manta/1.2.0; configManta.py --referenceFasta {params.genome} --runDir {output.dir}"+var
      print(cmd)
      shell(cmd)
      cmd="module load manta/1.2.0; module load python/2.7; {output.dir}/runWorkflow.py -m local -j {threads} -g 12"
      print(cmd)
      shell(cmd)