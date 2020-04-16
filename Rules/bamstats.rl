rule bamstats:
    input:  bam="{x}.dedup.bam",
            bed=config['project']['workpath']+"/exome_targets.bed"
    output: "{x}.dedup.bam.onTarget.bam_stats",
            temp("{x}.dedup.bam.onTarget.bam"),
            "{x}.dedup.bam.bam_stats"
    params: regions="exome_targets.bed",rname="pl:bamstats"
    shell:  "module load bamtools; module load samtools; module load bedtools; module load perl/5.18.4; perl Scripts/cal.on.target.pl --target_bed {params.regions} {input.bam}"