rule snpeff:
	input: "all.snp.vcf","all.indel.vcf"
        output: "all.snp.snpeff.vcf","all.indel.snpeff.vcf","stats_summary.snp.html","stats_summary.indel.html","all.snp.first.vcf","all.indel.first.vcf","all.snp.filter.vcf","all.indel.filter.vcf"
        params: snpeff=config['bin'][pfamily]['SNPEFF'],genome=config['references'][pfamily]['SNPEFF_GENOME'],snpsift=config['bin'][pfamily]['SNPSIFT'],effconfig=config['references'][pfamily]['SNPEFF_CONFIG'],firstanno=config['bin'][pfamily]['FIRSTANNO'],filtersift=config['references'][pfamily]['FILTERSIFT'],hapmap=config['references'][pfamily]['HAPMAP_DBNSFP'], regulomeDB=config['references'][pfamily]['REGULOMEDB'],rname="pl:snpeff"
        run:
            for mode in ["snp","indel"]:
                ifile="all.{}.vcf".format(mode)
                ohtml="stats_summary.{}.html".format(mode)
                ofirst="all.{}.first.vcf".format(mode)
                ofilter="all.{}.filter.vcf".format(mode)
                osnpeff="all.{}.snpeff.vcf".format(mode)
	        shell("module load snpEff;{params.snpeff} -v {params.genome} -c {params.effconfig} %s -interval {params.hapmap} interval {params.regulomeDB} -stats %s > %s; cat %s | python {params.firstanno}> %s; cat %s | {params.snpsift} filter \"{params.filtersift}\" > %s"%(ifile,ohtml,osnpeff,osnpeff,ofirst,ofirst,ofilter))	      
#	        shell("module load snpEff; {params.snpeff} -v {params.genome} -c {params.effconfig} {0} -stats {1} > {2}; cat {2} | python {params.firstanno}> {3}; cat {3} | {params.snpsift} filter \"{params.filtersift}\" > {4}".format(ifile,ohtml,osnpeff,ofirst,ofilter))

