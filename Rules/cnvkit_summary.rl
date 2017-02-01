rule cnvkit_summary:
    input: cnvkit=expand(config['project']['workpath']+"/cnvkit_out/{p}_calls.cns",p=pairs)
    output: heatmap=config['project']['workpath']+"/cnvkit_out/CNVkit_summary_heatmap.pdf"
    params: rname="pl:cnvkit_summary"
    shell: "module load cnvkit/0.8; cnvkit.py heatmap -d -o {output.heatmap} cnvkit_out/*.cns"