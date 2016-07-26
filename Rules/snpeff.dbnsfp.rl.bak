rule snpeff_dbnsfp:
	input: 'all.snp.filter.vcf','all.indel.filter.vcf'
	output: 'all.snp.dbnsfp.vcf','all.indel.dbnsfp.vcf'
	params: genome=config['references']['SNPEFF_GENOME'], snpsift=config['bin']['SNPSIFT'], effconfig=config['references']['SNPEFF_CONFIG'],selectDB=config['references']['SELECTDB'], db29=config['references']['DB29'],rname="pl:snpeff"
        run:
          for mode in ["snp","indel"]:
              ifile="all.{}.filter.vcf".format(mode)
              ofile="all.{}.dbnsfp.vcf".format(mode)
	      shell("module load snpEff; {params.snpsift} dbnsfp -f {params.selectDB} -v -db {params.db29} %s > %s"%(ifile,ofile))

