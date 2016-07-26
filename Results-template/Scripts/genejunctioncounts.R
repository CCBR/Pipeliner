args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
FILES <- args[2]
MINCOUNT <- as.numeric(args[4])
MINSAMPLES <- as.numeric(args[5])

GTFFILE <- args[3]
system(paste('grep -w gene ',GTFFILE,' | cut -f1,4,5,7,9 | sed -e "s/;/\\t/g" | cut -f1,2,3,4,5,7,9 | sed -e \'s/gene_id "//g\' -e \'s/gene_type "//g\' -e \'s/gene_name "//g\' | sed -e \'s/"//g\' > gencode.gtf.pre.bed ',sep=""))
pregtf=read.delim('gencode.gtf.pre.bed',header=F,stringsAsFactors=F)
filler1=rep(".",length(pregtf[,1]))
filler2=rep(".",length(pregtf[,1]))
id=paste(pregtf[,5],"-",pregtf[,7],sep="")
gtf=cbind(pregtf[,1],pregtf[,2],pregtf[,3],filler1,filler2,pregtf[,4],gsub("[[:blank:]]", "",id),gsub("[[:blank:]]", "",pregtf[,6]))
write.table(gtf,file="gencode.gtf.filler.bed",sep="\t",col.names=F,row.names=F,quote=F)


setwd(DIR)
myfiles=as.character(unlist(strsplit(FILES, split=" ")))

res=read.delim(myfiles[1],header=F,stringsAsFactors=F)
stranded=res[,4]
stranded=gsub("1","+",stranded)
stranded=gsub("2","-",stranded)
filler1 = rep(".",length(res[,1]))
filler2 = rep(".",length(res[,1]))
result=cbind(res[,1],res[,2],res[,3],filler1,filler2,stranded,res[,7]);
tmpoutfilename1=paste(substr(myfiles[1], 0, nchar(myfiles[1])-3),"bed.recoded.tab",sep="")
write.table(result,file=tmpoutfilename1,sep="\t",row.names=F,col.names=F,quote=F)
tmpoutfilename2=paste(substr(myfiles[1], 0, nchar(myfiles[1])-14),"_gencodeintersects.txt",sep="")
system(paste('module load bedtools/2.19.1; bedtools intersect -a ',"gencode.gtf.filler.bed",' -b ',tmpoutfilename1,' -wao -loj -s > ',tmpoutfilename2,sep=""))
tmpoutfilename3=paste(substr(myfiles[1], 0, nchar(myfiles[1])-14),"_genecounts.txt",sep="")
system(paste('awk \'{a[$7]+=$15}END{for(i in a){print i, a[i], OFS=\"\\t\"}}\' ',tmpoutfilename2," > ",tmpoutfilename3,sep=""))
tmp=read.delim(tmpoutfilename3,header=F,stringsAsFactors=F)
merged=cbind(id=tmp[,1],as.numeric(tmp[,2]))
colnames(merged)[2]=myfiles[1]
cat(tmpoutfilename1,"")
cat("1"," .. done \n")
for(i in seq(2, length(myfiles), by = 1))
{
	res=read.delim(myfiles[i],header=F,stringsAsFactors=F)
	stranded=res[,4]
	stranded=gsub("1","+",stranded)
	stranded=gsub("2","-",stranded)
	filler1 = rep(".",length(res[,1]))
	filler2 = rep(".",length(res[,1]))
	result=cbind(res[,1],res[,2],res[,3],filler1,filler2,stranded,res[,7]);
	tmpoutfilename1=paste(substr(myfiles[i], 0, nchar(myfiles[i])-3),"bed.recoded.tab",sep="")
	write.table(result,file=tmpoutfilename1,sep="\t",row.names=F,col.names=F,quote=F)
	tmpoutfilename2=paste(substr(myfiles[i], 0, nchar(myfiles[i])-14),"_gencodeintersects.txt",sep="")
	system(paste('module load bedtools/2.19.1; bedtools intersect -a ',"gencode.gtf.filler.bed",' -b ',tmpoutfilename1,' -wao -loj -s > ',tmpoutfilename2,sep=""))
	tmpoutfilename3=paste(substr(myfiles[i], 0, nchar(myfiles[i])-14),"_genecounts.txt",sep="")
	system(paste('awk \'{a[$7]+=$15}END{for(i in a){print i, a[i], OFS=\"\\t\"}}\' ',tmpoutfilename2," > ",tmpoutfilename3,sep=""))
	tmp=read.delim(tmpoutfilename3,header=F,stringsAsFactors=F)
	tomerge=cbind(id=tmp[,1],as.numeric(tmp[,2]))
	colnames(tomerge)[2]=myfiles[i]
	merged=merge(merged,tomerge,all=T)
	cat(tmpoutfilename1,"")
	cat(i," .. done \n")
}

for(i in seq(2, length(myfiles)+1, by = 1))
{
  ind=which(is.na(merged[,i]))
  if(!(length(ind)==0)){
  	merged[ind,i]=0
  }
  cat(i," .. done \n")
}

write.table(as.data.frame(merged),file="RawCountFile_genejunctions.txt",sep="\t",row.names=F,quote=F)
tofilter=read.delim("RawCountFile_genejunctions.txt",stringsAsFactors=F,row.names=1)


# Remove rows which do not pass filter
cat(MINCOUNT," ", MINSAMPLES, "checking..\n")
tot=colSums(tofilter)
MINCOUNT=(MINCOUNT/max(tot))*1e6
cat(MINCOUNT," ", MINSAMPLES, "checking..\n")
library("edgeR")
filter <- apply(cpm(tofilter), 1, function(x) length(x[x>MINCOUNT])>=MINSAMPLES)
#filter <- apply(myzero, 1, function(y) length(y[y>0])>=1)
filtered=tofilter[filter,]
write.table(as.data.frame(filtered),file="RawCountFile_genejunctions_filtered.txt",sep="\t",col.names=NA)
dim(filtered)

#Rscript genejunctioncounts.R '/scratch/gowdanb/rnatest/' 'SRR950078.p2.SJ.out.tab SRR950079.p2.SJ.out.tab SRR950080.p2.SJ.out.tab SRR950081.p2.SJ.out.tab SRR950082.p2.SJ.out.tab SRR950083.p2.SJ.out.tab SRR950084.p2.SJ.out.tab SRR950085.p2.SJ.out.tab SRR950086.p2.SJ.out.tab SRR950087.p2.SJ.out.tab' '/fdb/GENCODE/Gencode_human/release_19/gencode.v19.annotation.gtf' '5' '1'
