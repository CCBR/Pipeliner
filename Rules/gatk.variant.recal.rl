rule gatk_variant_recal:
	input: vcf="combined.gvcf"
        output: "all.snp.recal", "all.snp.tranches", "all.snp.plotting.R","all.indel.recal", "all.indel.tranches", "all.indel.plotting.R"
        params: GATK=config['bin'][pfamily]['GATK'], REF=config['references'][pfamily]['GENOME'],rname="pl:varrecal"
        run: 
          recal={"snp":"-mode SNP -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -resource:hapmap,known=false,training=true,truth=true,prior=15.0 " + config['references'][pfamily]['HAPMAP'] + " -resource:omni,known=false,training=true,truth=true,prior=12.0 " + config['references'][pfamily]['OMNI'] + " -resource:1000G,known=false,training=true,truth=false,prior=10.0 " + config['references'][pfamily]['B1K'] + " -resource:dbsnp,known=true,training=false,truth=false,prior=2.0  " + config['references'][pfamily]['SNP138'],"indel":"-mode INDEL -an DP -an FS -an MQRankSum -an ReadPosRankSum -resource:mills,known=true,training=true,truth=true,prior=12.0 " + config['references'][pfamily]['INDELSITES']}
          for mode in recal.keys():
              ofile="all."+mode+".recal"
              tfile="all."+mode+".tranches"
              rfile="all."+mode+".plotting.R"
              shell("{params.GATK} -T VariantRecalibrator -R {params.REF} -input {input.vcf} "+recal[mode]+" -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile "+ofile+" -tranchesFile "+tfile+" -rscriptFile "+rfile)

