suppressMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-d", "--degout", type="character", required=TRUE,
                    help="Reformmated DEG out file from limma/edgeR/DESeq2")
parser$add_argument("-c", "--gene2coord", type="character", required=TRUE,
                    help="Gene to coordinate file")
parser$add_argument("-g", "--genome", type="character", required=TRUE,
                    help="Genome .. either hg19/hg38/mm9/mm10/Mmul8.0.1/canFam3/hg38_30")
parser$add_argument("-f", "--fdr", type="double", default=0.05,
                    help="FDR cutoff to use")
parser$add_argument("-o", "--outfileprefix", type="character", required=TRUE,
                    help="DEG type output file prefix: limma/edgeR/DESeq2")

args <- parser$parse_args()


for (f in c(args$degout,args$gene2coord))
if (! file.exists(f)) {
  stop(paste("File does not exist:",f))
}

if (! args$genome %in% c("hg19","hg38","hg38_30","mm9","mm10","Mmul8.0.1","canFam3")) {
  stop("Only hg19/hg38/mm9/mm10/Mmul8.0.1/canFam3/hg38_30 genomes are supported!")
}


suppressMessages(library("karyoploteR"))
suppressMessages(library("BSgenome.Mmusculus.UCSC.mm9"))
suppressMessages(library("BSgenome.Mmusculus.UCSC.mm10"))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg19"))
suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38"))
suppressMessages(library("BSgenome.Mmulatta.UCSC.rheMac8"))
suppressMessages(library("BSgenome.Cfamiliaris.UCSC.canFam3.masked"))

if (args$genome=="Mmul8.0.1"){args$genome="rheMac8"}

deseqout=read.delim(args$degout)
dim(deseqout)
fdr_filter=deseqout$fdr < args$fdr
positive_lfc_filter=deseqout$log2fc>1
negative_lfc_filter=deseqout$log2fc < -1
table(( negative_lfc_filter | positive_lfc_filter ) & fdr_filter )
deseqout_filtered=deseqout[( ( negative_lfc_filter | positive_lfc_filter ) & fdr_filter ),]

coordinates=read.delim(args$gene2coord)
deseqout_filtered_w_coord=merge(deseqout_filtered,coordinates,by.x="gene",by.y="gene_name")
dim(deseqout_filtered_w_coord)

if(nrow(deseqout_filtered_w_coord)==0){
  print("Warning: No significant differentially expressed genes found! Please try increasing the FDR cutoff.")
}

genome=args$genome
chrs=c()
maxchrs=0
if (genome %in% c("hg19","hg38","hg38_30")) {maxchrs=22}
if (genome %in% c("mm10","mm9")) {maxchrs=19}
if (genome %in% c("rheMac8")) {maxchrs=20}
if (genome %in% c("canFam3")) {maxchrs=38}

for (i in seq(1,maxchrs)) {chrs=c(chrs,paste("chr",i,sep=""))}
chrs=c(chrs,"chrX")
if (! genome %in% c("canFam3")){chrs=c(chrs,"chrY")}

y=round(length(chrs)/2)
a=chrs[seq(1,y)]
b=chrs[seq(y+1,length(chrs))]
chrs_subsets=list(a,b)

deseqout_filtered_w_coord=deseqout_filtered_w_coord[deseqout_filtered_w_coord$chr %in% chrs,]
dim(deseqout_filtered_w_coord)

scale_limit=max(abs(floor(fivenum(deseqout_filtered_w_coord$log2fc)[2])),ceiling(fivenum(deseqout_filtered_w_coord$log2fc)[4]))+1
neg_scale_limit = -1 * scale_limit

pos_strand_up_triangle=deseqout_filtered_w_coord[ (deseqout_filtered_w_coord$strand=="+" & deseqout_filtered_w_coord$log2fc>scale_limit),]
pos_strand_up_dot=deseqout_filtered_w_coord[ (deseqout_filtered_w_coord$strand=="+" & deseqout_filtered_w_coord$log2fc>0 & deseqout_filtered_w_coord$log2fc<scale_limit),]

neg_strand_up_triangle=deseqout_filtered_w_coord[ (deseqout_filtered_w_coord$strand=="-" & deseqout_filtered_w_coord$log2fc>scale_limit),]
neg_strand_up_dot=deseqout_filtered_w_coord[ (deseqout_filtered_w_coord$strand=="-" & deseqout_filtered_w_coord$log2fc>0 & deseqout_filtered_w_coord$log2f<scale_limit),]

pos_strand_down_triangle=deseqout_filtered_w_coord[ (deseqout_filtered_w_coord$strand=="+" & deseqout_filtered_w_coord$log2fc<neg_scale_limit),]
pos_strand_down_dot=deseqout_filtered_w_coord[ (deseqout_filtered_w_coord$strand=="+" & deseqout_filtered_w_coord$log2fc>neg_scale_limit & deseqout_filtered_w_coord$log2fc<0),]

neg_strand_down_triangle=deseqout_filtered_w_coord[ (deseqout_filtered_w_coord$strand=="-" & deseqout_filtered_w_coord$log2fc<neg_scale_limit),]
neg_strand_down_dot=deseqout_filtered_w_coord[ (deseqout_filtered_w_coord$strand=="-" & deseqout_filtered_w_coord$log2fc>neg_scale_limit & deseqout_filtered_w_coord$log2fc<0),]



if (nrow(pos_strand_up_triangle)>0) {pos_strand_up_triangle$log2fc=scale_limit}
if (nrow(neg_strand_up_triangle)>0) {neg_strand_up_triangle$log2fc=scale_limit}
if (nrow(pos_strand_down_triangle)>0) {pos_strand_down_triangle$log2fc=neg_scale_limit}
if (nrow(neg_strand_down_triangle)>0) {neg_strand_down_triangle$log2fc=neg_scale_limit}

for (i in seq(1,length(chrs_subsets))) {
  chrs2=unlist(chrs_subsets[i])
  png(paste(args$outfileprefix,"_karyoplot",i,".png",sep=""), width = 10, height = 6, units = 'in', res = 1600)
  kp <- plotKaryotype(genome=genome, plot.type=2, chromosomes = chrs2)

  kpDataBackground(kp, data.panel = 1, r0=0, r1=0.9)
  kpDataBackground(kp, data.panel = 2, r0=0, r1=0.9)
  kpAxis(kp, ymin=1, ymax=scale_limit, r0=0, r1=1, col="gray50", cex=0.25, data.panel = 1)
  kpAxis(kp, ymin=-1, ymax=neg_scale_limit, r0=0, r1=1, col="gray50", cex=0.25, data.panel = 2)
  kpPoints(kp, chr=pos_strand_up_dot$chr, x=pos_strand_up_dot$coord, y=pos_strand_up_dot$log2fc, col="red", pch=".", cex=0.5,data.panel = 1, ymin=1, ymax=scale_limit)
  kpPoints(kp, chr=pos_strand_down_dot$chr, x=pos_strand_down_dot$coord, y=pos_strand_down_dot$log2fc, col="darkgreen", pch=".", cex=0.5,data.panel = 2, ymin=-1, ymax=neg_scale_limit)
  kpPoints(kp, chr=neg_strand_up_dot$chr, x=neg_strand_up_dot$coord, y=neg_strand_up_dot$log2fc, col="red", pch=".", cex=0.5,data.panel = 1, ymin=1, ymax=scale_limit)
  kpPoints(kp, chr=neg_strand_down_dot$chr, x=neg_strand_down_dot$coord, y=neg_strand_down_dot$log2fc, col="darkgreen", pch=".", cex=0.5,data.panel = 2, ymin=-1, ymax=neg_scale_limit)

  kpPoints(kp, chr=pos_strand_up_triangle$chr, x=pos_strand_up_triangle$coord, y=pos_strand_up_triangle$log2fc, col="red", pch=2, cex=0.25,data.panel = 1, ymin=1, ymax=scale_limit)
  kpPoints(kp, chr=pos_strand_down_triangle$chr, x=pos_strand_down_triangle$coord, y=pos_strand_down_triangle$log2fc, col="darkgreen", pch=6, cex=0.25,data.panel = 2, ymin=-1, ymax=neg_scale_limit)
  kpPoints(kp, chr=neg_strand_up_triangle$chr, x=neg_strand_up_triangle$coord, y=neg_strand_up_triangle$log2fc, col="red", pch=2, cex=0.25,data.panel = 1, ymin=1, ymax=scale_limit)
  kpPoints(kp, chr=neg_strand_down_triangle$chr, x=neg_strand_down_triangle$coord, y=neg_strand_down_triangle$log2fc, col="darkgreen", pch=6, cex=0.25,data.panel = 2, ymin=-1, ymax=neg_scale_limit)
  dev.off()
}




