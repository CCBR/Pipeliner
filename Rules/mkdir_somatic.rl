rule mkdir_somatic:
        input:  expand("{s}.recal.bam", s=samples),expand("{s}.recal.bam.bai", s=samples)
        output: temp("QC/decoy"),
                mutect2_dir=config['project']['workpath']+"/mutect2_out",
                strelka_dir=config['project']['workpath']+"/strelka_out",
                mutect_dir=config['project']['workpath']+"/mutect_out",
        params: rname="pl:mkdir"
        shell:  "echo \'decoy\' > {output}; mkdir -p strelka_out; mkdir -p strelka_out/oncotator_out; mkdir -p strelka_out/mutsigCV_out; mkdir -p cnvkit_out; mkdir -p mutect_out; mkdir -p mutect_out/oncotator_out; mkdir -p mutect_out/mutsigCV_out; mkdir -p theta_out; mkdir -p conpair_out; mkdir -p germline_vcfs; mkdir -p mutect2_out; mkdir -p mutect2_out/oncotator_out; mkdir -p mutect2_out/mutsigCV_out; mkdir -p mutect2_out/chrom_files; mkdir -p manta_out"