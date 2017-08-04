library("maftools")
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
FILE1 <- args[2]
FILE2 <- args[3]
FILE3 <- args[4]
setwd(DIR)
mymaf <- read.maf(FILE1)
plotmafSummary(mymaf,file=FILE2)
pdf(FILE3)
oncoplot(mymaf,writeMatrix=TRUE,showTumorSampleBarcodes=TRUE)
dev.off()