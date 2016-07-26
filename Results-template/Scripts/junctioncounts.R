args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
FILES <- args[2]
MINCOUNT <- as.numeric(args[3])
MINSAMPLES <- as.numeric(args[4])

setwd(DIR)
myfiles=as.character(unlist(strsplit(FILES, split=" ")))

res=read.delim(myfiles[1],header=F)
id=paste(res[,1],"-",res[,2],"-",res[,3],"-",res[,4],sep="")
result=cbind(id,res[,7]); colnames(result)[2]=myfiles[1]
cat("1"," .. done \n")
# read files
 for(i in seq(2, length(myfiles), by = 1))
#  for(i in seq(2, 2, by = 1))
{
  res=read.delim(myfiles[i],header=F)
  id=paste(res[,1],"-",res[,2],"-",res[,3],"-",res[,4],sep="")
  temp=cbind(id,res[,7]); colnames(temp)[2]=myfiles[i]
  result=merge(result,temp,all=T)
  cat(i," .. done \n")
 }

write.table(result,file="RawCountFile_junctions.txt",sep="\t",row.names=F)
myresult=read.delim("RawCountFile_junctions.txt",row.names=1)

for(i in seq(1, length(myfiles), by = 1))
{
  ind=which(is.na(myresult[,i]))
  if(!(length(ind)==0)){
  	myresult[ind,i]=0
  }
  cat(i," .. done \n")
}

myzero=myresult
write.table(myzero,file="RawCountFileZero_junctions.txt",sep="\t",col.names=NA)

# Remove rows which do not pass filter
cat(MINCOUNT," ", MINSAMPLES, "checking..\n")
tot=colSums(myzero)
MINCOUNT=(MINCOUNT/max(tot))*1e6
cat(MINCOUNT," ", MINSAMPLES, "checking..\n")
library("edgeR")
filter <- apply(cpm(myzero), 1, function(x) length(x[x>MINCOUNT])>=MINSAMPLES)
#filter <- apply(myzero, 1, function(y) length(y[y>0])>=1)
myzero=myzero[filter,]
write.table(as.data.frame(myzero),file="RawCountFile_junctions_filtered.txt",sep="\t",col.names=NA)
dim(myzero)