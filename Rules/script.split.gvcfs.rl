rule script_split_gvcfs:
    input: "combined.gvcf"
    output: expand("all.{m}.vcf",m=["snp","indel"])
    params: rname="pl:splitgvcfs"
    shell: "Scripts/split.combined.py"
