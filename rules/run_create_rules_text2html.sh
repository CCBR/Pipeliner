RULES_DIR=/Users/mcintoshc/Desktop/ExomePipeline-master/Rules

cd /Users/mcintoshc/Documents/REPOSITORIES/CEM/ExomePipeline-gh-pages/rules
for RULE in "add.sample.tag" "all-bam2recal" "all-custom.combined.template" "all-custom" "all-custom.template" "all-exomeseq-germline-partial" "all-exomeseq-germline-realign" "all-exomeseq-germline" "all-exomeseq-pairs" "all-exomeseq-somatic" "all-initialqc" "all-wgslow" "annovar" "annovar2" "avia.make.bed.somatic" "avia" "bicseq" "bwa.indexref" "bwa.mem" "bwa.pe" "fastqc.fastq" "fastqc.trimmed" "gatk.apply.variant.recal" "gatk.apply.vqsr" "gatk.combine.gvcfs" "gatk.genotype.gvcfs" "gatk.genotyping" "gatk.haplotype.caller.recal" "gatk.haplotype.caller" "gatk.merge.somatic.vcfs" "gatk.mutect2" "gatk.realign.2" "gatk.realign" "gatk.recal" "gatk.variant.annnotator" "gatk.variant.recal" "make.somatic.network" "mutect" "ngsqc" "novo.sort" "novocraft.novoalign" "novocraft.sort" "picard.headers.old" "picard.headers" "picard.headers2" "picard.index.bam" "picard.markdups" "qualimap" "rules.db" "samtools.flagstats.dedup" "samtools.flagstats" "samtools.mpileup" "samtools.sam2bam" "script.batchgvcf" "script.bzip2" "script.checkqc" "script.coverage.qc" "script.fastqc.summary" "script.split.gvcfs" "snpeff_somatic" "snpeff.dbnsfp" "snpeff.new" "snpeff" "somatic.germline.calls" "star.align.2" "star.align" "star.generate" "svdetect" "trimmomatic" "varscan"
do
	python create_rules_text2html.py \
	--rule_file $RULES_DIR/$RULE.rl \
	--rule_name $RULE \
	--out_file $RULE.html
done