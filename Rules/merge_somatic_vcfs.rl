rule merge_somatic_vcfs:
    input: mutect=config['project']['workpath']+"/mutect_out/{x}.FINAL.vcf",
           strelka=config['project']['workpath']+"/strelka_out/{x}_FINAL.vcf",
           mutect2=config['project']['workpath']+"/mutect2_out/{x}.FINALmutect2.vcf",
    output: mergedvcf=config['project']['workpath']+"/merged_somatic_callers/merged_somatic.vcf",
    params: regions="exome_targets.bed",gres="lscratch:100",gatk=config['bin'][pfamily]['GATK'],genome=config['references'][pfamily]['GENOME'],snpsites=config['references'][pfamily]['SNPSITES'],rname="CombineVariants"
    shell: "module load GATK/3.8-0; java -Xmx48g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID -jar $GATK_JAR -T CombineVariants -R {params.genome} -nt 8 --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED --genotypemergeoption PRIORITIZE --rod_priority_list mutect2,mutect,strelka -o {output.mergedvcf} --variant:mutect {input.mutect} --variant:strelka {input.strelka} --variant:mutect2 {input.mutect2}"