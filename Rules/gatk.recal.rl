rule gatk_recal:
      input:  "{x}.realign.bam"
      output: bam="{x}.recal.bam",
              re=temp("{x}_recal_data.bam.grp"),
              re2=temp("{x}_recal_data2.bam.grp")
#              plots="{X}.plots.bam.pdf"
      params: gatk=config['bin'][pfamily]['GATK'],
              genome=config['references'][pfamily]['GENOME'],
              indelsites=config['references'][pfamily]['INDELSITES'],
              snpsites=config['references'][pfamily]['SNPSITES'],rname="pl:recal"
      shell:  """
              {params.gatk} -T BaseRecalibrator -I {input} -R {params.genome} -knownSites {params.snpsites} -knownSites {params.indelsites} -o {output.re}; 
              {params.gatk} -T PrintReads -R {params.genome} -I {input} -BQSR {output.re} -o {output.bam}; 
              {params.gatk} -T BaseRecalibrator -I {output.bam} -R {params.genome} -knownSites {params.snpsites} -knownSites {params.indelsites} -BQSR {output.re} -o {output.re2}
              """
# -plots {output.plots}
