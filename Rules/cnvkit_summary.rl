rule cnvkit_summary:
    input: cnvkit=expand("cnvkit_out/{p}_calls.cns",p=samples)
    output: heatmap="cnvkit_out/CNVkit_summary_heatmap.pdf",
#            metrics="cnvkit_out/CNVkit_summary.metrics"
    params: tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],rname="pl:cnvkit_summary"
    run:
      fl=os.popen("ls cnvkit_out/*.cns").read().split()
      var=" "+" ".join(fl)
      cmd="module load cnvkit/0.8; cnvkit.py heatmap -d -o {output.heatmap} "+var
      print(cmd)
      shell(cmd)