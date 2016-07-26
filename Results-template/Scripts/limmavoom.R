## library('edgeR')
library('statmod')
library('RColorBrewer') 
library('gplots')
library('reshape') 
library('ggplot2')
library('limma')
library('geneplotter')

## grab args
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
FILE1 <- args[2]
FILE2 <- args[3]
CONTRASTS <- args[4]

setwd(DIR)
## Sys.setenv("DISPLAY"=":0.0")
## options(device=NULL)
# read files
sampleinfo=read.delim(FILE1)
x = read.delim(FILE2,row.names=1)
# sampleFiles=as.character(sampleinfo[,2])
Group <- factor(sampleinfo$condition)
design=model.matrix(~0+Group)
contras=unlist(strsplit(CONTRASTS, split=" "))        
cat(contras,"\t",length(contras),"\t",contras[1],"\t",contras[2],"\n",file="readcontraLimma.txt")
cons=c()
for(i in seq(1, length(contras), by = 2))
{{
cons=c(cons,paste(contras[i],"-",contras[i+1],sep=""))
}}
#print(x)
#print(design)
png("VoomPlot.png")
v1 <- voom(as.matrix(x),design,plot=TRUE,normalize="quantile")
dev.off()
sf = v1$E/log2((x/colSums(x))*1000000)
write.table(sf,file="LimmaVoom_scaling_factors.txt",row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
png("HistLimmavoomNormFilter.png")
df.n <- melt(as.data.frame(v1$E))
print(ggplot(df.n) + geom_density(aes(x = value,colour = variable)) + labs(x = NULL) + theme(legend.position='top'))
dev.off()
## MDS
png("Limma_MDS.png")
## MDS <- plotMDS(v1,xlim=c(-5,5),ylim=c(-5,5),cex=1,pch=20)
plotMDS(v1,xlim=c(-5,5),ylim=c(-5,5),cex=1,pch=20)
## png("Limma_MDS.png")
shortname=paste(substr(colnames(v1$E),1,22))
## text(MDS, labels=shortname, cex=0.5, pos=1)
## print(MDS)
## dev.copy(png, paste("Limma_MDS.png"))
dev.off()
## DEG
nb=length(contras)/2
colnames(design) <- levels(Group)
fit <- lmFit(v1,design)
contrast.matrix <- makeContrasts(contrasts=cons,levels=design)
fitb <- contrasts.fit(fit, contrast.matrix)
ebayes.fit=eBayes(fitb)
## 
for (i in 1:nb)
{{
all.genes.con = topTable(ebayes.fit, coef = i, number=nrow(ebayes.fit))
## generate Volcano plot
jpeg(paste("Limma_",cons[i],"_volcano.jpeg",sep=""),quality=100) 
plot(all.genes.con$logFC,-log10(all.genes.con$adj.P.Val),cex=0.1,xlab="Log Fold-Change",ylab="-log10 Adj P-Value",main=paste('Volcano Plot for ',cons[i],sep=""))
t=which(all.genes.con$adj.P.Val<0.05 & abs(all.genes.con$logFC)>=1 )
points(all.genes.con$logFC[t],-log10(all.genes.con$adj.P.Val[t]),col="red",pch=20,cex=0.5)
dev.off()
#MAplot <- plot(ebayes.fit,coef=i)
#print(MAplot)
#dev.copy(png, paste(cons[i],"_MAplot_Limma_old.png",sep=""))
#dev.off()
dataf=data.frame("m"=all.genes.con$AveExpr,"fc"=all.genes.con$logFC,"sig"=all.genes.con$adj.P.Val<0.05)
png(paste(cons[i],"_MAplot_Limma_v2.png",sep=""))
plotMA(dataf,log="",main=cons[i],ylim=range(all.genes.con$logFC))
dev.off()
all.genes.con$FC <- ifelse(all.genes.con$logFC<0, -1/(2^all.genes.con$logFC), 2^all.genes.con$logFC)
write.table(all.genes.con,file=paste("Limma_deg_",cons[i],"_all_genes.txt",sep=""),sep="\t",col.names=NA)
}}
#

