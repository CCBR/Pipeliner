rule mkdir_somatic:
        input:  expand("{s}.recal.bam", s=samples)
        output: mutect2_dir=config['project']['workpath']+"/mutect2_out",
                mutect_dir=config['project']['workpath']+"/mutect_out",
                vardict_dir=config['project']['workpath']+"/vardict_out",
                merged_dir=config['project']['workpath']+"/merged_somatic_variants",
                mut2chrom_dir=config['project']['workpath']+"/mutect2_out/chrom_files",
                manta_dir=config['project']['workpath']+"/manta_out"
        params: rname="pl:mkdir"
        shell:  "mkdir -p germline_vcfs; mkdir -p mutect2_out; mkdir -p mutect2_out/oncotator_out; mkdir -p mutect2_out/mutsigCV_out; mkdir -p mutect2_out/chrom_files; mkdir -p mutect_out; mkdir -p mutect_out/oncotator_out; mkdir -p mutect_out/mutsigCV_out; mkdir -p vardict_out; mkdir -p vardict_out/oncotator_out; mkdir -p vardict_out/mutsigCV_out; mkdir -p merged_somatic_variants; mkdir -p merged_somatic_variants/oncotator_out; mkdir -p merged_somatic_variants/mutsigCV_out; mkdir -p manta_out"