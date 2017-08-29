#!/usr/local/apps/R/gcc_4.9.1/3.2.3/bin/Rscript
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

#read first file for gene length data
myfirstfiles = gsub("count", "count.info", myfiles[1])
res=read.delim(myfirstfiles,header=T, comment.char = '#')
colnames(res)[7]=as.character(myfiles[1]) 

#read all other files and merge to the first one
for(i in seq(2, length(myfiles), by = 1))
{{
  temp=read.delim(myfiles[i],header=T)
  colnames(temp)[2]=as.character(myfiles[i]) 
  res=merge(res,temp)
}}

#read annotation flie and merge all together
gene_name=read.delim(ANNOTATE,header=F,sep=" ")
res2=merge(gene_name,res,by.x=1,by.y=1)

#reformat and output raw counts
res2=cbind(Ensembl_id=res2[,1],Gene_symbol=res2[,3],res2[,-c(1:9)])
write.table(as.data.frame(res2[,-c(3)]),file="RawCount_genes_unfiltered.txt",sep="\t",row.names=F) 

#using DEGlist 
mydata<-DGEList(counts=res2[,-c(1:3)],genes=res2[,c("Ensembl_id","Gene_symbol","Length")])

#apply filter to counts
MINCOUNT<-5
MINSAMPLES<-2
val1=as.numeric(MINCOUNT)
val2=as.numeric(MINSAMPLES)

tot=colSums(res2[,-c(1:3)])
val1=(val1/max(tot))*1e6

filter <- apply(cpm(mydata), 1, function(x) length(x[x>val1])>=val2)
res=mydata[filter,,keep.lib.sizes=FALSE]

#output filtered raw counts, still output old format for the sake of pipeline, will change it soon.

write.table(data.frame(cbind(symbol=paste(res$genes[,1],"|",res$genes[,2],sep=""),res$counts)),file="RawCountFile_genes_filtered.txt",sep="\t", row.names = F)

#draw before png
png("HistBeforenormFilter.png")
df.m <- melt(as.data.frame(res$counts))
print(ggplot(df.m) + geom_density(aes(x = value, colour = variable)) + labs(x = NULL) + theme(legend.position='top') + scale_x_log10())
dev.off() 

#calculate CPM_TMM and output
tmm_y <- calcNormFactors(res,method="TMM")
ndata= cpm(tmm_y,log=FALSE,normalized.lib.sizes=TRUE)
write.table(data.frame(res$genes[,-c(3)], ndata),file="CPM_TMM_counts.txt",sep="\t",row.names=F)

#calculate RPKM_TMM and output
ndata = rpkm(tmm_y,res$genes$Length, log=FALSE,normalized.lib.sizes=TRUE)
write.table(data.frame(res$genes[,-c(3)], ndata),file="RPKM_TMM_counts.txt",sep="\t",row.names=F)

#calculate unfiltered CPM_TMM and output
tmm_y <- calcNormFactors(mydata,method="TMM")
ndata= cpm(tmm_y,log=FALSE,normalized.lib.sizes=TRUE)
write.table(data.frame(mydata$genes[,-c(3)], ndata),file="CPM_TMM_unfiltered_counts.txt",sep="\t",row.names=F)

#calculate unfiltered RPKM_TMM and output
ndata = rpkm(tmm_y,res$genes$Length, log=FALSE,normalized.lib.sizes=TRUE)
write.table(data.frame(mydata$genes[,-c(3)], ndata),file="RPKM_TMM_unfiltered_counts.txt",sep="\t",row.names=F)
