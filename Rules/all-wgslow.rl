rule all_wgslow:
    input: "combined.vcf",
           "full_annot.txt.zip",
           "variants.database",
           "sample_network.bmp",
           "combined.snpeff.vcf"
    output: 
    params: rname="final"
    shell:  """
             mv *.out slurmfiles/

            """