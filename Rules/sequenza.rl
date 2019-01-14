rule sequenza:
    input: freeccnvs="freec_out/pass1/{x}.recal.bam_CNVs",
    output: segments="sequenza_out/{x}/{x}"+"_segments.txt",
            fit="sequenza_out/{x}/{x}"+"_alternative_solutions.txt",
            seqz=temp("sequenza_out/{x}/{x}.seqz.gz"),
            binseqz=temp("sequenza_out/{x}/{x}.bin100.seqz.gz"),
    params: dir=config['project']['workpath']+"/sequenza_out/{x}",tumorsample=lambda wildcards: config['project']['pairs'][wildcards.x][1],normalsample=lambda wildcards: config['project']['pairs'][wildcards.x][0],gc=config['references'][pfamily]['SEQUENZAGC'],rname="pl:sequenza"
    threads: 8
    shell: "mkdir -p sequenza_out; mkdir -p sequenza_out/{params.tumorsample}; module load sequenza-utils/2.2.0; module load samtools/1.9; gzip -c freec_out/pass1/{params.tumorsample}/{params.normalsample}.recal.bam_minipileup.pileup > sequenza_out/{params.tumorsample}/{params.normalsample}.recal.bam_minipileup.pileup.gz; gzip -c freec_out/pass1/{params.tumorsample}/{params.tumorsample}.recal.bam_minipileup.pileup > sequenza_out/{params.tumorsample}/{params.tumorsample}.recal.bam_minipileup.pileup.gz; sequenza-utils bam2seqz -p -gc {params.gc} -n sequenza_out/{params.tumorsample}/{params.normalsample}.recal.bam_minipileup.pileup.gz -t sequenza_out/{params.tumorsample}/{params.tumorsample}.recal.bam_minipileup.pileup.gz | gzip > {output.seqz}; sequenza-utils seqz_binning -w 100 -s {output.seqz} | gzip > {output.binseqz}; module load R/3.5; Rscript Scripts/run_sequenza.R {output.binseqz} {params.dir} {threads} {params.normalsample}+{params.tumorsample}"