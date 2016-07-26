library('DESeq2')
library('RColorBrewer') 
library('gplots')
library('reshape') 
library('ggplot2')

## grab args
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
FILE1 <- args[2]
FILE2 <- args[3]
ANNOTATE <- args[4]
CONTRASTS <- args[5]

#
setwd(DIR)
## Sys.setenv("DISPLAY"=":0.0")
sampleinfo=read.delim(FILE1)
sampleFiles=as.character(sampleinfo[,2])
x = read.delim(FILE2,row.names=1)
## read annotation file
## ann=read.delim(ANNOTATE)
#
ddsHTSeq<-DESeqDataSetFromMatrix(countData=x,colData=sampleinfo, design=~condition)
dds<-DESeq(ddsHTSeq)
ndata=as.data.frame(counts(dds,normalized=TRUE))
colnames(ndata)=colnames(x)
write.table(ndata,file="Deseq2_normalized_counts.txt",sep="\t",col.names=NA)
png("HistDesq2normFilter.png")
df.m <- melt(as.data.frame(ndata))
print(ggplot(df.m) + geom_density(aes(x = value, colour = variable)) + labs(x = NULL) + theme(legend.position='top') + scale_x_log10())
dev.off() 
#
contras=unlist(strsplit(CONTRASTS, split=" "))        
cat(contras,"\t",length(contras),"\t",contras[1],"\t",contras[2],"\n",file="readcontra.txt")
for(i in seq(1, length(contras), by = 2))
{{
res<-results(dds,contrast=c("condition",as.character(contras[i]),as.character(contras[i+1])))
res<-res[order(res$padj),]
res1=as.data.frame(res)
restmp=res1
restmp$FoldChange <- ifelse(restmp$log2FoldChange<0, -1/(2^restmp$log2FoldChange), 2^restmp$log2FoldChange)
write.table(restmp,file=paste("DEG_",contras[i],"_vs_",contras[i+1],".txt",sep=""),sep="\t",col.names=NA) 
x=res1$log2FoldChange[which(!is.na(res1$log2FoldChange))] 
png(paste("MAplot_",contras[i],"_vs_",contras[i+1],".png",sep=""))
plotMA(res,ylim=range(x),main=paste("MAplot_",contras[i],"_vs_",contras[i+1],sep=""))
##dev.copy(png,paste("MAplot_",contras[i],"_vs_",contras[i+1],".png",sep=""))
dev.off()
}}
## transformation
rld <- rlogTransformation(dds, blind=TRUE)
rldm=assay(rld)
colnames(rldm)=colnames(x)
write.table(rldm,file="Deseq2_normalized_rld.txt",sep="\t",col.names=NA)
## clustering
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),paste(condition,sampleFiles , sep=" : "))
#if you just want the conditions use this line : rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
png("deseq2_heatmaps_samplebysample.png")
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(16, 16))
## dev.copy(png,"deseq2_heatmaps_samplebysample.png")
dev.off()
## plotMA(dds,ylim=c(-2,2),main="DESeq2 MAplot")
## dev.copy(png,"deseq2_MAplot.png")
## dev.off()
## pca
png("deseq2_pca.png")
print(plotPCA(rld, intgroup=c("condition")))
## dev.copy(png,"deseq2_pca.png")
dev.off()
png("deseq2_pca_details.png")
print(plotPCA(rld, intgroup=c("condition","fileName")))
## dev.copy(png,"deseq2_pca_details.png")
dev.off()
