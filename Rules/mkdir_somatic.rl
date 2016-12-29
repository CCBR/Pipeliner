rule mkdir_somatic:
        input:  expand("{s}.dedup.bam", s=samples)
        output: mutect_dir="mutect_out",
                cnvkit_dir="cnvkit_out",
                strelka_dir="strelka_out",
                delly_dir="delly_out",
                theta_dir="theta_out",
                conpair_dir="conpair_out"
        params: rname="pl:mkdir"
        shell:  "mkdir -p strelka_out; mkdir -p strelka_out/oncotator_out; mkdir -p strelka_out/mutsigCV_out; mkdir -p cnvkit_out; mkdir -p delly_out; mkdir -p mutect_out; mkdir -p mutect_out/oncotator_out; mkdir -p mutect_out/mutsigCV_out; mkdir -p theta_out; mkdir -p conpair_out; mkdir -p germline_vcfs; mkdir -p mutect2_out; mkdir -p mutect2_out/chrom_files"