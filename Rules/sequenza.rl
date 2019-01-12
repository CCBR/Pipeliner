rule sequenza:
    input: normalpileup="freec_out/pass1/{x}/"+lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam_minipileup.pileup",
           tumorpileup="freec_out/pass1/{x}/"+lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam_minipileup.pileup",
    output: normalpileupgz=temp("sequenza_out/{x}/"+lambda wildcards: config['project']['pairs'][wildcards.x][0]+".recal.bam_minipileup.pileup.gz"),
            tumorpileupgz=temp("sequenza_out/{x}/"+lambda wildcards: config['project']['pairs'][wildcards.x][1]+".recal.bam_minipileup.pileup.gz"),
            segments="sequenza_out/{x}/{x}"+"_segments.txt",
            fit="sequenza_out/{x}/{x}"+"_alternative_solutions.txt",
            seqz=temp("sequenza_out/{x}/"+lambda wildcards: config['project']['pairs'][wildcards.x][1]+".seqz.gz"),
            binseqz=temp("sequenza_out/{x}/"+lambda wildcards: config['project']['pairs'][wildcards.x][1]+".bin100.seqz.gz"),
    params: pair={x},dir=config['project']['workpath']/sequenza_out/{x},tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],rname="pl:sequenza"
    threads: 8
    shell: "mkdir -p sequenza_out
            mkdir -p sequenza_out/{params.pair}
            module load sequenza-utils/2.2.0
            module load samtools/1.9
            gzip -c {input.normalpileup} > {output.normalpileupgz}
            gzip -c {input.tumorpileup} > {output.tumorpileupgz}
            sequenza-utils bam2seqz -p -gc {params.gc} -n {output.normalpileupgz} -t {output.tumorpileupgz} | gzip > {output.seqz}
            sequenza-utils seqz_binning -w 100 -s {output.seqz} | gzip > {output.binseqz}
            module load R/3.5
            Rscript Scripts/run_sequenza.R {output.binseqz} {params.dir} {threads} {params.pair}"