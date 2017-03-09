rule snpeff_somatic:
        input: mutect=config['project']['workpath']+"/mutect_out/merged_somatic.vcf",
               strelka=config['project']['workpath']+"/strelka_out/merged_somatic.vcf"
        output: mutect=config['project']['workpath']+"/mutect_out/merged_somatic_snpEff.vcf",
                strelka=config['project']['workpath']+"/strelka_out/merged_somatic_snpEff.vcf"
        params: snpeff=config['bin'][pfamily]['SNPEFF'],regions="exome_targets.bed",genome=config['references'][pfamily]['SNPEFF_GENOME'],rname="pl:snpeff_somatic"
        shell: "module load snpEff; {params.snpeff} -v {params.genome} -interval {params.regions} -cancer -canon -csvStats snpeff.mutect.stats.csv -stats snpeff.mutect.stats.html -cancerSamples pairs {input.mutect} > {output.mutect}; {params.snpeff} -v {params.genome} -interval {params.regions} -cancer -canon -csvStats snpeff.strelka.stats.csv -stats snpeff.strelka.stats.html -cancerSamples pairs {input.strelka} > {output.strelka}"