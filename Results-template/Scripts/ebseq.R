args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
FILES <- args[2]
SAMTAB <- args[3]
CONTRASTS <- args[4]
RSEMREF <- args[5]
RSEM <- args[6]

#DIR = '/scratch/dwheeler/test'
#FILES = '12284_N_S38.rsem.isoforms.results 12284_T_S37.rsem.isoforms.results 21238_N_S42.rsem.isoforms.results 21238_T_S41.rsem.isoforms.results 23709_N_S40.rsem.isoforms.results 23709_T_S39.rsem.isoforms.results'
#SAMTAB = 'sampletable.txt'
#CONTRASTS = 'N T'
#RSEMREF = "/data/CCBR/dev/RNA-Seq-pipeline/ver3/rsemref123/human_gencode"
#RSEM = "/data/CCBR/apps/rsem1.2.23/RSEM-1.2.23/"

setwd(DIR)
myfiles=as.character(unlist(strsplit(FILES, split=" ")))

#library("EBSeq", lib.loc=PIPERLIB)

s2c<-read.table(file.path(paste(DIR,"/",SAMTAB,sep="")), header=TRUE, stringsAsFactors = FALSE)
s2c<-dplyr::select(s2c, sample=sampleName, condition)
s2c$sample<-gsub("\\..*","",s2c$sample)

contras=unlist(strsplit(CONTRASTS, split=" "))

for(f in seq(1, length(contras), by = 2)){

	s2c_pair = s2c[which(s2c$condition==as.character(contras[f])|s2c$condition==as.character(contras[f+1])),]
	pairname = paste(contras[f],"_vs_",contras[f+1],sep="")
	cat(paste(pairname,'\n',sep=''))
	files = myfiles[which(gsub("\\..*","",myfiles) %in% s2c_pair$sample)]
	fileorder = vector()
	c1 = 0
	for(fl in files){
		if(s2c_pair$condition[match(gsub("\\..*","",fl),s2c_pair$sample)] == contras[f]){
			fileorder = c(fileorder, fl)
			c1 = c1+1

		}
	}
	c2 = 0
	for(fl in files){
		if(s2c_pair$condition[match(gsub("\\..*","",fl),s2c_pair$sample)] == contras[f+1]){
			fileorder = c(fileorder, fl)
			c2 = c2 + 1

		}
	}
	isofiles = paste(fileorder,collapse=" ")

	system(paste(RSEM,"/rsem-generate-data-matrix ",isofiles," > ",pairname,".isoform.counts.matrix",sep=""))

	system(paste(RSEM,"/rsem-run-ebseq --ngvector ",RSEMREF,".transcripts.ngvec ",pairname,".isoform.counts.matrix ",c1,",",c2," ",pairname,".isoform.EBSeq",sep=""))

}

write.csv(s2c, file="ebseq_completed.txt")
