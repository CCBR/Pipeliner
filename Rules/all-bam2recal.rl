rule all_bam2recal:
    input:  expand("{s}"+".recal.bam",s=samples),
    output:
