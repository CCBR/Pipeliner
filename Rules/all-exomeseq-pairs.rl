rule all_exomeseq_pairs:
    input:  expand("{s}"+".pileup.bam",s=samples),
#            expand("{s}.coverage.done",s=samples),
            expand("{s}.snp",s=pairs),
            expand("{s}.snp.vcf.av",s=pairs),
            expand("{s}.uniqcov",s=samples),
#            expand("QC/{s}.qualimapReport.pdf",s=samples),
#            expand("QC/{s}.flags",s=samples),
#            expand("QC/{s}.dedup.flags",s=samples)
    output: 
#    params: s=expand("{s}",s=samples)
#    shell: "Rscript QC/align_qc.R {params.s}"
