library('reshape') 
library('ggplot2')
library('edgeR')
library('DESeq2')
library('tidyverse')

writegzfile <- function(m,f) {
  m=as.data.frame(m)
  m$id=rownames(m)
  m=separate(data=m,col=id,into=c('ensID','geneName'),sep="\\|",remove=TRUE)
  m=m %>% select('ensID','geneName',everything())
  write.table(m,file=gzfile(f),sep="\t",row.names = FALSE,quote=F)
} 

args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
FILES <- args[2]
ANNOTATE <- args[3]
SAMPLETABLE <- args[4]

#SAMPLETABLE="fullsampletable.txt"
#DIR="~/Desktop/Temp/ccbr842/RNASeq/RSEM_filtering"
#FILES=c("iv_4T1C1_1.RSEM.genes.results iv_4T1C1_2.RSEM.genes.results iv_NF1_1.RSEM.genes.results iv_NF1_2.RSEM.genes.results iv_TSC1_1.RSEM.genes.results iv_TSC1_2.RSEM.genes.results iv_Tgfbr_2_1.RSEM.genes.results iv_Tgfbr_2_2.RSEM.genes.results")
#ANNOTATE="annotate.genes.txt"

MINCOUNT=0.5
MINSAMPLES=0.5


setwd(DIR)
x=read.table(SAMPLETABLE,header = T,sep="\t")
myfiles=as.character(unlist(strsplit(FILES, split=" ")))
res=read.delim(myfiles[1],header=T)[,c(1,5)] 
colnames(res)[2]=as.character(myfiles[1]) 
for(i in seq(2, length(myfiles), by = 1))
{{
temp=read.delim(myfiles[i],header=T)[,c(1,5)]
colnames(temp)[2]=as.character(myfiles[i]) 
res=merge(res,temp)
}}

gene_name=read.delim(ANNOTATE,header=F,sep=" ")
res2=merge(gene_name,res,by.x=1,by.y=1)
res3=cbind(symbol=paste(res2[,1],"|",res2[,3],sep=""),res2[,-c(1,2,3,4,5)])
write.table(as.data.frame(res3),file="RawCountFile_RSEM_genes.txt",sep="\t",row.names=F,quote = F) 
mydata=read.delim("RawCountFile_RSEM_genes.txt",row.names=1)
colnames(mydata)=gsub('\\..*$','',colnames(mydata))
mydata=ceiling(mydata)

writegzfile(cpm(mydata),"RSEM_CPM_counts.txt.gz")

groups=levels(x$condition)
G1=groups[1]
g1_samples=(x$condition==G1)
ng1=max(1,floor(length(g1_samples[g1_samples==TRUE])*MINSAMPLES))
CPM_CUTOFF=MINCOUNT
mydata1=mydata[,g1_samples]
k_g1=rowSums(cpm(mydata1)>CPM_CUTOFF)>=ng1
k=k_g1
table(k)

for(i in seq(2,length(levels(x$condition)))){
  Gi=groups[i]
  gi_samples=(x$condition==Gi)
  ngi=max(1,floor(length(gi_samples[gi_samples==TRUE])*MINSAMPLES))
  mydatai=mydata[,gi_samples]
  k_gi=rowSums(cpm(mydatai)>CPM_CUTOFF)>=ngi
  k=k|k_gi
  print(table(k))
}

res=mydata[k,]
res2=res
res2$symbol=rownames(res2)
res2=res2 %>% select('symbol',everything())
write.table(res2,file="RawCountFile_RSEM_genes_filtered.txt",row.names = F,quote = F,sep="\t")
# png("RSEM_HistBeforenormFilter.png")
# df.m <- melt(as.data.frame(res))
# print(ggplot(df.m) + geom_density(aes(x = value, colour = variable)) + labs(x = NULL) + theme(legend.position='top') + scale_x_log10())
# dev.off() 
y = DGEList(counts=res)
## Normalization TMM ------------------------------------------------------------
## method = =c("TMM","RLE","upperquartile","none")
y <- calcNormFactors(y,method="TMM")
ndata= cpm(y,log=FALSE,normalized.lib.sizes=TRUE)
## save it 
writegzfile(ndata,"RSEM_CPM_TMM_counts.txt.gz")
	## unfiltered normalization
y2 = DGEList(counts=mydata)
y2 <- calcNormFactors(y2,method="TMM")
ndata2= cpm(y2,log=FALSE,normalized.lib.sizes=TRUE)
## save it 
writegzfile(ndata2,"RSEM_CPM_TMM_unfiltered_counts.txt.gz")
rlogres=rlog(as.matrix(res),blind=TRUE)
rownames(rlogres)=rownames(res)
writegzfile(rlogres,"RSEM_rlog_counts.txt.gz")
