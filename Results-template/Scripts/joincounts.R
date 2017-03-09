# -------------------------------------------------------------------------------------
# joining the counts from Subread with Overlap
#
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
FILES <- args[2]
FILES2 <- args[3]
STARSTRANDCOL <- args[4]
setwd(DIR)
myfiles=as.character(unlist(strsplit(FILES, split=" ")))
res=read.delim(myfiles[1],header=T)
colnames(res)[1]="gene"
colnames(res)[2]=as.character(myfiles[1]) 
# remove the last 5 statistics lines ... 
# nr=dim(res)[1]
# res=res[-c((nr-4):nr),]
#
for(i in seq(2, length(myfiles), by = 1))
{{
temp=read.delim(myfiles[i],header=T)
colnames(temp)[1]="gene"
colnames(temp)[2]=as.character(myfiles[i]) 
res=merge(res,temp)
}}
write.table(as.data.frame(res),file="RawCountFileOverlap.txt",sep="\t",row.names=F) 
# -------------------------------------------------------------------------------------
# joining the counts from star
# here we assume unstranded so take column TWO => should a PARAMETER in json file 
#
# vcol=2
vcol=as.numeric(STARSTRANDCOL)
myfiles2=as.character(unlist(strsplit(FILES2, split=" ")))
res2=read.delim(myfiles2[1],header=F)[,c(1,vcol)]
colnames(res2)[1]="gene"
colnames(res2)[2]=as.character(myfiles2[1]) 
# remove the first  4  statistics lines ... 
# 
res2=res2[-c(1:4),]
#
for(i in seq(2, length(myfiles2), by = 1))
{{
temp=read.delim(myfiles2[i],header=F)[,c(1,vcol)]
temp=temp[-c(1:4),]
colnames(temp)[1]="gene"
colnames(temp)[2]=as.character(myfiles2[i]) 
res2=merge(res2,temp)
}}
write.table(as.data.frame(res2),file="RawCountFileStar.txt",sep="\t",row.names=F) 
