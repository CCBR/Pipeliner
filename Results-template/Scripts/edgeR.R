library('edgeR')
library('statmod')
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
# read files
sampleinfo=read.delim(FILE1)
x = read.delim(FILE2,row.names=1)
# sampleFiles=as.character(sampleinfo[,2])
## read annotation file
## ann=read.delim(ANNOTATE)
# DGElist object --------------------------------------------------------------
condition = as.factor(sampleinfo$condition)
y = DGEList(counts=x,group=condition)
## Normalization TMM ------------------------------------------------------------
## method = =c("TMM","RLE","upperquartile","none")
y <- calcNormFactors(y,method="TMM")
# y$samples
png("libdistrib.png")
barplot(y$samples$lib.size*1e-6,main="Library size distribution", names= strsplit(colnames(y$counts),".star.count.txt"), ylab="Library size (millions)",las=2,cex.names=0.8)
dev.off()
## MDS plots ----------------------------------------------------------------------
# both pairewise (leading)
png("MDS_bcv.png")
# print(y)
plotMDS(y, method="bcv", , main="MDS plot bcv")
dev.off()
png("MDS_logFC.png")
plotMDS(y, method="logFC" , main="MDS plot logFC") ## plotMDS(y) default
dev.off()
# plotMDS(y, method="logFC",gene.selection="common", main="MDS plot common")
## estimating common and tagwise dispersions -----------------------------------------
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y) #default trend: moveingave
## plotting
png("BCVplot.png")
plotBCV(y,main="BCV plot")
dev.off()
## differentially expressed genes ---------------------------------------------------
contras=unlist(strsplit(CONTRASTS, split=" "))        
cat(contras,"\t",length(contras),"\t",contras[1],"\t",contras[2],"\n",file="readcontra-edgeR.txt")
for(i in seq(1, length(contras), by = 2))
{{
deg<-exactTest(y,c(as.character(contras[i+1]),as.character(contras[i])))
# 
n=dim(y$counts)[1]
tt=topTags(deg, n=n)
res1 = as.data.frame(tt)
#
## res1=cbind(Ensembl.Gene.ID=substr(rownames(res1),1,18),id.ver=rownames(res1),res1)
## final= merge(res1,ann,all.x=TRUE)
##final=final[order(final$FDR),]
final=res1[order(res1$FDR),]
final$FC <- ifelse(final$logFC<0, -1/(2^final$logFC), 2^final$logFC)
write.table(final,file=paste("DEG_EdgeR_",contras[i],"_vs_",contras[i+1],".txt",sep=""),sep="\t",col.names=NA)
#  like MAplot
deg1sel <- decideTestsDGE(deg, p=0.05, adjust="BH")
detags <- rownames(y)[as.logical(deg1sel)]
png(paste("Smearplot_",contras[i],"_vs_",contras[i+1],".png",sep=""))
plotSmear(deg, de.tags=detags,main= paste("Smearplot FDR<0.05 ",contras[i],"_vs_",contras[i+1],sep=""))
abline(h = c(-2, 2), col = "blue")
dev.off()
# 
}}
## transformation
ylog2=cpm(y,log=TRUE,normalized.lib.sizes=TRUE,prior.count=2) # prior count like avelogcpm
ndata= cpm(y,log=FALSE,normalized.lib.sizes=TRUE)*1e6
## save it
write.table(ylog2,file="edgeR_normalized_counts_log.txt",sep="\t",col.names=NA) 
write.table(ndata,file="edgeR_normalized_counts.txt",sep="\t",col.names=NA)
png("HistEdgeRnormFilter.png")
df.m <- melt(as.data.frame(ndata))
print(ggplot(df.m) + geom_density(aes(x = value, colour = variable)) + labs(x = NULL) + theme(legend.position='top') + scale_x_log10())
dev.off()         
## clustering / heatmap
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
distylog2=dist(t(ylog2))
mat = as.matrix(distylog2)
# rownames(mat) <- colnames(mat)
png("edgeR_heatmaps_samplebysample.png")
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(16, 16))
## dev.copy(png,"edgeR_heatmaps_samplebysample.png")
dev.off()
#pca
pr2=prcomp(t(ylog2))
png("edgeR_prcomp.png")
# biplot(pr2)
plot(pr2$x[,1],pr2$x[,2],col="red", main="PCA plot using prcomp and Logcpm data")
text(pr2$x[,1],pr2$x[,2], labels=colnames(ylog2), cex=0.7, pos=4)
dev.off()
