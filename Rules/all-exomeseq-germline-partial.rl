rule all_exomeseq_germline_partial:
    input:  "combined.gvcf",
            expand("all.{type}.dbnsfp.vcf", type=["snp","indel"])    
#            expand("all.{m}.recal",m=["snp","indel"]),
#            expand("QC/{s}.qualimapReport.pdf",s=samples),
#            expand("QC/{s}.flags",s=samples),
#            expand("QC/{s}.dedup.flags",s=samples)
    output: 
