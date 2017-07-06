read.table("temp.database", sep="\t", stringsAsFactors=F, check.names=F, header=T) -> vcf
vcf <- vcf[,!duplicated(colnames(vcf))]
write.table(vcf, "temp.database", quote=F, row.names=F, col.names=T, sep="\t")
