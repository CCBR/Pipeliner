library('reshape') 
library('ggplot2')
library('edgeR')
library('limma')
library('tidyverse')
library('DESeq2')

writegzfile <- function(m,f) {
  m=as.data.frame(m)
  m$id=rownames(m)
  m=separate(data=m,col=id,into=c('ensID','geneName'),sep="\\|",remove=TRUE)
  m=m %>% select('ensID','geneName',everything())
  write.table(m,file=gzfile(f),sep="\t",row.names = FALSE,quote=F)
}

args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
FILE1 <- args[2]
MINCOUNT <- args[3]
MINSAMPLES <- args[4]
ANNOTATE <- args[5]
FILE2 <- args[6]

setwd(DIR)

myfiles=as.character(unlist(strsplit(FILE1, split=" ")))

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
# res2=cbind(Ensembl_id=res2[,1],Gene_symbol=res2[,3],res2[,-c(1:9)])
res3=cbind(symbol=paste(res2[,1],"|",res2[,3],sep=""),res2[,-c(1:9)])
colnames(res3)=gsub('\\..*$','',colnames(res3))
colnames(res3)=gsub('.*/','',colnames(res3))
write.table(as.data.frame(res3[,-c(2)]),file="RawCountFile_Subread_genes.txt",sep="\t",row.names=F,quote = F)

#using DEGlist 
res4=res3[,-c(2)]
rownames(res4)=res4$symbol
res4$symbol=NULL
# mydata<-DGEList(counts=res3[,-c(1:2)],genes=res3[,c("symbol","Length")])
mydata<-DGEList(res4)

#apply filter to counts
#MINCOUNT<-5
#MINSAMPLES<-2
val1=as.numeric(MINCOUNT)
val2=as.numeric(MINSAMPLES)

filter <- apply(cpm(mydata), 1, function(x) length(x[x>val1])>=val2)
res=mydata[filter,,keep.lib.sizes=FALSE]

#output filtered raw counts, still output old format for the sake of pipeline, will change it soon.
res2=as.data.frame(res$counts)
res2$symbol=rownames(res2)
res2=res2 %>% select('symbol',everything())
write.table(res2,file="RawCountFile_Subread_genes_filtered.txt",row.names = F,quote = F,sep="\t")

#draw before png
png("HistBeforenormFilter.png")
df.m <- melt(as.data.frame(res$counts))
print(ggplot(df.m) + geom_density(aes(x = value, colour = variable)) + labs(x = NULL) + theme(legend.position='top') + scale_x_log10())
dev.off() 

#calculate CPM_TMM and output
tmm_y <- calcNormFactors(res,method="TMM")
ndata= cpm(tmm_y,log=FALSE,normalized.lib.sizes=TRUE)
writegzfile(ndata,"Subread_CPM_TMM_counts.txt.gz")

rlogres=rlog(as.matrix(res),blind=TRUE)
rownames(rlogres)=rownames(res)
writegzfile(rlogres,"Subread_rlog_counts.txt.gz")
