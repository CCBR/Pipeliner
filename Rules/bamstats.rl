rule bamstats:
    input:  "{x}.dedup.bam"
    output: "{x}.dedup.bam.onTarget.bam_stats",
            "{x}.dedup.bam.onTarget.bam",
            "{x}.dedup.bam.bam_stats"
    params: regions=config['references'][pfamily]['REFFLAT'],rname="pl:bamstats"
    shell:  "module load bamtools; module load samtools; module load bedtools; perl Scripts/cal.on.target.pl --target_bed {params.regions} {input}"