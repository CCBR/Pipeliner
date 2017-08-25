rule mkdir_somatic:
        input:  expand("{s}.recal.bam", s=samples)
        output: mutect2_dir=config['project']['workpath']+"/mutect2_out",
                mut2chrom_dir=config['project']['workpath']+"/mutect2_out/chrom_files",
                cnvkit_dir=config['project']['workpath']+"/cnvkit_out",
                delly_dir=config['project']['workpath']+"/delly_out",
                manta_dir=config['project']['workpath']+"/manta_out"
        params: rname="pl:mkdir"
        shell:  "mkdir -p cnvkit_out; mkdir -p delly_out; mkdir -p germline_vcfs; mkdir -p mutect2_out; mkdir -p mutect2_out/oncotator_out; mkdir -p mutect2_out/mutsigCV_out; mkdir -p mutect2_out/chrom_files; mkdir -p manta_out"