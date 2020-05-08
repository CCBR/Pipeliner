rule generate_QC_table:
    input: expand("{s}.dedup.bam.onTarget.bam_stats",s=samples), expand("{s}.dedup.bam.bam_stats",s=samples), expand("QC/{s}_run_trimmomatic.err",s=samples), expand("QC/{s}.qualimapReport/genome_results.txt",s=samples), expand("{s}.sorted.txt",s=samples), "multiqc_report.html"
    output: config['project']['id']+"_"+config['project']['flowcellid']+".xlsx"
    params: project=config['project']['id'],flowcell=config['project']['flowcellid'],dir=config['project']['workpath'],rname="pl:QC_table"
    shell: "module load perl/5.18.4; perl Scripts/CollectPipelineStats2Tab_v2.3.pl -p {params.project} -f {params.flowcell} -d {params.dir} -r 3 -e 2; perl Scripts/Tab2Excel_v2.3.pl -i {params.project}_{params.flowcell} -r 3"