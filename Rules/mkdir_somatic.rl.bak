rule mkdir_somatic:
        input:  
        output: mutect2="mutect2",
                mutect="mutect",
                cnvkit="cnvkit",
                strelka="strelka",
                delly="delly"
        params: rname="pl:mkdir"
        shell:  "mkdir {output.delly}; mkdir {output.mutect2}; mkdir {output.mutect}; mkdir mutect/oncotator; mkdir mutect/mutsigCV; mkdir {output.strelka}; mkdir {output.cnvkit}; mkdir strelka/oncotator; mkdir strelka/mutsigCV; mkdir mutect2/oncotator; mkdir mutect2/mutsigCV"
