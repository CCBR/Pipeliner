rule fusion_summary:
    input: expand("starfusion/oncofuse/{x}/{x}.oncofuse.output", x=samples),
           expand("fusioncatcher/oncofuse/{x}/{x}.oncofuse.output"),
    output: dir="summary",
            file="summary/summaryCounts.tsv"
    params: rname="pl:summary"
    shell: "Scripts/fusionSummary.sh"