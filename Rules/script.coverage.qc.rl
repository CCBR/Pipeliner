rule script_coverage_qc:
    input: "{x}.pileup.bam"
    output: hist = "{x}.coverage.done"
    params: rname="pl:coverage"
    shell:  "cat {input} |cut -f1,2,4 > out && Rscript QC/coverage.R {wildcards.x};touch {wildcards.x}.coverage.done"

