rule mkdir_somatic:
        input:  expand("{s}.dedup.bam", s=samples)
        output: mutect_dir="mutect_out",
                cnvkit_dir="cnvkit_out",
                strelka_dir="strelka_out",
                delly_dir="delly_out"
        params: rname="pl:mkdir"
        shell:  "mkdir strelka_out cnvkit_out delly_out mutect_out mutect_out/oncotator_out mutect_out/mutsigCV_out strelka_out/oncotator_out strelka_out/mutsigCV_out theta_out conpair_out germline_vcfs"