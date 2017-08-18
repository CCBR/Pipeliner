rule cnvkit_summary:
    input: cnvkit=expand(config['project']['workpath']+"/cnvkit_out/{s}_calls.cns",s=samples)
    output: heatmap=config['project']['workpath']+"/cnvkit_out/CNVkit_summary_heatmap.pdf"
    params: rname="pl:cnvkit_summary"
    shell: "module load cnvkit/0.8; cnvkit.py heatmap -d -o {output.heatmap} cnvkit_out/*.cns; module load R/3.4.0_gcc-6.2.0; Rscript Scripts/cnplot.R cnvkit_out"