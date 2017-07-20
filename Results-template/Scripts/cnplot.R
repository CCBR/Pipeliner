#adding an argument
args <- commandArgs(TRUE)
args <- setwd(args[1])
#Open directory and looping through the files with .cns extention  

dataset<-NULL
file_list <- list.files(pattern ="*.cns")


for (file in file_list){
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t")
  }
    if (exists("dataset")){
    temp_dataset <-read.table(file, header=TRUE, sep="\t")
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  
}


n<-length(file_list)


#Select data and columns needed 
data <- dataset[,c("chromosome","start","end","cn","cn1","cn2")]
data$chromosome<-toupper(data$chromosome)
data$chromosome <- gsub('CHR', '', data$chromosome)
data <- data[data$chromosome != 'X' & data$chromosome != 'Y',]
data$chromosome <- as.numeric(as.character(data$chromosome))
data<-data[order(data$chromosome, data$start ),]


#Shifting the positions to take account the chromosomes 
posmin <- tapply(data$end ,data$chromosome, min);
posmax <- tapply(data$end,data$chromosome, max);
posshift <- head(c(0,cumsum(as.numeric(posmax))),-1);
names(posshift) <- levels(data$chromosome)
genpos <- data$start + posshift[data$chromosome];

#Using the shifts 
data$newStart <- genpos
data$newEnd <- data$newStart + abs(data$start - data$end)
data$points <-mapply(seq, data$newStart, data$newEnd,50000)

#Seperating data based on type 
gain<-data[data$cn >2 & data$cn <5,]
loss<-data[data$cn ==1,]
loh<-data[(data$cn1 == 0 | data$cn2 == 0) & data$cn >0,]
loh<-loh[!is.na(loh$cn1) & !is.na(loh$cn2),]
homozygous<-data[data$cn == 0,]
hla<-data[data$cn >4,]
breaks<-length(seq(min(unlist(data$points)),max(unlist(data$points)),50000))

#Starting the plot
pdf(file="CNplot.pdf",width = 16,height = 8,compress = TRUE)

h1<-hist(unlist(gain$points),breaks = breaks,plot = FALSE)
h2<-hist(unlist(loss$points),breaks= breaks, plot = FALSE)
h3<-hist(unlist(loh$points),breaks= breaks, plot = FALSE)
if (nrow(homozygous) >0 ) {
  h4<-hist(unlist(homozygous$points),breaks= breaks, plot = FALSE)
}

if (nrow(hla)>0) {
  h5<-hist(unlist(hla$points),breaks= breaks, plot = FALSE)
}


h1$counts = (h1$counts/n)*100
h2$counts = -(h2$counts/n)*100
h3$counts = -(h3$counts/n)*100
if (nrow(homozygous) >0 ) {
  h4$counts = -(h4$counts/n)*250
}

if (nrow(hla) >0 ) {
  h5$counts =(h5$counts/n)*250
                }
begend <-c(posshift,max(data$newEnd))
chromN<-length(begend)-1

#Y=c(h1$counts,h2$counts,h3$counts)
#ymax = max(Y)
#ymin = min(Y)


plot(h1, ylim=c(-100, 100), col="indianred1",yaxt='n',xaxt='n',lty="blank",xlab = NULL,main = "",ylab = "")
lines(h3, col="gray",lty="blank")
lines(h2, col="lightblue2",lty="blank")
if (nrow(homozygous) >0 ) {
  lines(h4, col="blue4",lty="blank")
  }
if (nrow(hla) >0 ) {
  lines(h5, col="red4",lty="blank")
}
axis(1,begend,line = 1,col="black",col.ticks="black",col.axis="black",labels = c(seq(1:chromN),''))
axis(2,c(-80,-60,-40,-20,0,20,40,60,80),labels=c('80%','60%','40%','20%',0,'20%','40%','60%','80%'))
axis(4,c(-75,-50,-25,0,25,50,75),labels = c('30%','20%','10%',0,'10%','20%','30%'))
mtext(text = 'Gain', side = 2,col = 'indianred1',las=1, adj=0,line=3, at=90)
mtext(text = 'Loss', side = 2,col = 'lightblue3',las=1, adj=0,line=3, at=-95)
mtext(text = 'LOH', side = 2,col = 'darkgray',las=1, adj=0,line=3, at=-90)
mtext(text = 'HLA', side = 4,col = 'red4',las=1, adj=0,line=-1, at=85 )
mtext(text = 'Deletion', side = 4,col = 'blue3',las=1, adj=0,line=-1, at=-85)
abline(v=begend[1:22],col='black',lty=5)
graphics.off()




