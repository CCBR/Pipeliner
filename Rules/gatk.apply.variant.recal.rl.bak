rule apply_variant_recal:
	input: vcf="combined.gvcf","all.snp.recal", "all.snp.tranches","all.indel.recal", "all.indel.tranches"
        output: "all.vqsr.snp.vcf", "all.vqsr.indel.vcf" 
        params: GATK=config['bin']['GATK'], REF=config['references']['GENOME'],,rname="pl:recalvar"
        run: 
          for mode in ["snp","indel"]:
              rfile="all."+mode+".recal"
              tfile="all."+mode+".tranches"
              ofile="all.vqsr."+mode+".vcf"
	      shell ("{params.GATK} -T ApplyRecalibration -R {params.REF} -input {input.vcf} --ts_filter_level 99.0 -tranchesFile "+tfile+" -recalFile "+rfile+" -mode "+mode+" -o "+ofile)





