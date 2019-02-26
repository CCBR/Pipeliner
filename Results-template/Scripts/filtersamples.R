library(limma)
library(edgeR)
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
G1 <- args[2]
G2 <- args[3]
CPM_CUTOFF <- as.numeric(args[4])
MINSAMPLES <- as.numeric(args[5])
SAMPLETABLE <- args[6]
RAWCOUNTSTABLE <- args[7]
FILTEREDRAWCOUNTSTABLE <- args[8]

x=read.table(SAMPLETABLE,header = T,sep="\t")

samples=(x$condition==G1 | x$condition==G2)
x1=x[samples,]

g1_samples=(x1$condition==G1)
ng1=min(length(g1_samples[g1_samples==TRUE]),MINSAMPLES)
g2_samples=(x1$condition==G2)
ng2=min(length(g2_samples[g2_samples==TRUE]),MINSAMPLES)

y=read.table(RAWCOUNTSTABLE,header = T,sep = "\t")
y=y[,samples]
row.names(y)=y$symbol
y=y[,-which(names(y) %in% ("symbol"))]
y=ceiling(y)

y_g1=y[,g1_samples]
y_g2=y[,g2_samples]
# kcount_g1=rowSums(y_g1>10)>=ng1
# kcount_g2=rowSums(y_g2>10)>=ng2
k_g1=rowSums(cpm(y_g1)>CPM_CUTOFF)>=ng1
k_g2=rowSums(cpm(y_g2)>CPM_CUTOFF)>=ng2
table(k_g1)
table(k_g2)
# table(kcount_g1)
# table(kcount_g2)
k=k_g1|k_g2
table(k)
# k=k|kcount_g1|kcount_g2
# table(k)
y=y[k,]
y=cbind(symbol=row.names(y),y)
write.table(y,file=FILTEREDRAWCOUNTSTABLE,sep="\t",col.names=T,quote=F,row.names = F)

