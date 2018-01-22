library('reshape') 
library('ggplot2')
library('edgeR')

args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
FILES <- args[2]
MINCOUNT <- args[3]
MINSAMPLES <- args[4]
ANNOTATE <- args[5]

setwd(DIR)
myfiles=as.character(unlist(strsplit(FILES, split=" ")))
res=read.delim(myfiles[1],header=T)[,c(1,5)]
# colnames(res)[1]="gene"
colnames(res)[2]=as.character(myfiles[1]) 
# remove the last 5 statistics lines ... 
# nr=dim(res)[1]
# res=res[-c((nr-4):nr),]
#
for(i in seq(2, length(myfiles), by = 1))
{{
temp=read.delim(myfiles[i],header=T)[,c(1,5)]
#colnames(temp)[1]="gene"
colnames(temp)[2]=as.character(myfiles[i]) 
res=merge(res,temp)
}}
gene_name=read.delim(ANNOTATE,header=F,sep=" ")
res2=merge(gene_name,res,by.x=1,by.y=1)
res3=cbind(symbol=paste(res2[,1],"|",res2[,3],sep=""),res2[,-c(1,2,3,4,5)])
write.table(as.data.frame(res3),file="RawCountFile_RSEM_genes.txt",sep="\t",row.names=F) 
#
mydata=read.delim("RawCountFile_RSEM_genes.txt",row.names=1)
#  rounding
mydata=round(mydata)
val1=as.numeric(MINCOUNT)
val2=as.numeric(MINSAMPLES)
# cat(val1," ", val2, "checking..\n",file="check.txt")
## filter <- apply(mydata, 1, function(x) length(x[x>val1])>=val2)
## res=mydata[filter,]
#tot=colSums(mydata)
#val1=(val1/max(tot))*1e6
# cat(val1," ", val2, "checking..\n",file="check2.txt")
filter <- apply(cpm(mydata), 1, function(x) length(x[x>val1])>=val2)
res=mydata[filter,]
write.table(as.data.frame(res),file="RawCountFile_RSEM_genes_filtered.txt",sep="\t",col.names=NA)
png("RSEM_HistBeforenormFilter.png")
df.m <- melt(as.data.frame(res))
print(ggplot(df.m) + geom_density(aes(x = value, colour = variable)) + labs(x = NULL) + theme(legend.position='top') + scale_x_log10())
dev.off() 
y = DGEList(counts=res)
## Normalization TMM ------------------------------------------------------------
## method = =c("TMM","RLE","upperquartile","none")
y <- calcNormFactors(y,method="TMM")
ndata= cpm(y,log=FALSE,normalized.lib.sizes=TRUE)
## save it 
write.table(ndata,file="RSEM_CPM_TMM_counts.txt",sep="\t",col.names=NA)
	## unfiltered normalization
y2 = DGEList(counts=mydata)
y2 <- calcNormFactors(y2,method="TMM")
ndata2= cpm(y2,log=FALSE,normalized.lib.sizes=TRUE)
## save it 
write.table(ndata2,file="RSEM_CPM_TMM_unfiltered_counts.txt",sep="\t",col.names=NA)

