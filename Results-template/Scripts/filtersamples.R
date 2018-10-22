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
#write.table(x1,file=OUTSAMPLETABLE,sep="\t",col.names=T,quote=F,row.names = F)


y=read.table(RAWCOUNTSTABLE,header = T,sep = "\t")
row.names(y)=y$symbol
y=y[,-which(names(y) %in% ("symbol"))]
#y=y[,samples]
y=ceiling(y)
#colnames(y)=x1$label
k=rowSums(cpm(y)>CPM_CUTOFF)>=MINSAMPLES
# table(k)
y=y[k,]

y1=cbind(symbol=row.names(y),y)
write.table(y1,file=FILTEREDRAWCOUNTSTABLE,sep="\t",col.names=T,quote=F,row.names = F)
