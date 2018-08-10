#Example Usage: Rscript limma_edgeR_DESeq2_venn.R --contrast 'group1-group2' --limma limma_DEG_{group1}-{group2}_all_genes.txt --edgeR edgeR_DEG_{group1}-{group2}_all_genes.txt --DESeq2 DESeq2_deg_{group1}_vs_{group2}.txt

rm(list=ls())
library(Vennerable)
library(argparse)
library(dplyr)

filter_by_fc_fdr <-function(df,fdr_cutoff,upreg_fc_cutoff,downreg_fc_cutoff) {
  df=filter(df, fdr!='NA' & fc!='NA')
  keep1=df$fdr<=fdr_cutoff
  keep2=df$fc>=upreg_fc_cutoff
  keep3=df$fc<=downreg_fc_cutoff
  keep <- (keep1) & (keep2 | keep3)
  return(df[keep,])
}


parser <- ArgumentParser()

parser$add_argument("-c", "--contrast", type="character", required=TRUE,
                    help="contrast for differential expression")

parser$add_argument("-l", "--limma", type="character", required=TRUE,
                    help="reformatted limma output filename ")

parser$add_argument("-e", "--edgeR", type="character", required=TRUE,
                    help="reformatted edgeR output filename ")

parser$add_argument("-D", "--DESeq2", type="character", required=TRUE,
                    help="reformatted DESeq2 output filename ")

parser$add_argument("-f", "--fdr_cutoff", type="double", default=0.05,
                    help="FDR cutoff")

parser$add_argument("-u", "--upreg_fc_cutoff", type="double", default=1.5,
                    help="upregulated fold change cutoff")

parser$add_argument("-d", "--downreg_fc_cutoff", type="double", default=-1.5,
                    help="upregulated fold change cutoff")

args <- parser$parse_args()

limma=read.table(args$limma,header = TRUE)
edger=read.table(args$edgeR,header=TRUE)
deseq2=read.table(args$DESeq2,header=TRUE)

fdr_cutoff=args$fdr_cutoff
upreg_fc_cutoff=args$upreg_fc_cutoff
downreg_fc_cutoff=args$downreg_fc_cutoff


limma_filt=filter_by_fc_fdr(limma,fdr_cutoff,upreg_fc_cutoff,downreg_fc_cutoff)
edger_filt=filter_by_fc_fdr(edger,fdr_cutoff,upreg_fc_cutoff,downreg_fc_cutoff)
deseq2_filt=filter_by_fc_fdr(deseq2,fdr_cutoff,upreg_fc_cutoff,downreg_fc_cutoff)

dim(limma_filt)
limma_filt
dim(edger_filt)
edger_filt
dim(deseq2_filt)
deseq2_filt

x.genes=limma_filt$ensid_gene
y.genes=edger_filt$ensid_gene
z.genes=deseq2_filt$ensid_gene

vennD=Venn(SetNames = c("limma","edgeR","DESeq2"), 
           Weight=c(`100`=length(setdiff(x.genes,union(y.genes,z.genes))),
                    `010`=length(setdiff(y.genes,union(x.genes,z.genes))),
                    `110`=length(setdiff(intersect(x.genes,y.genes),z.genes)),
                    `001`=length(setdiff(z.genes,union(x.genes,y.genes))),
                    `101`=length(setdiff(intersect(x.genes,z.genes),y.genes)),
                    `011`=length(setdiff(intersect(z.genes,y.genes),x.genes)),
                    `111`=length(intersect(intersect(x.genes,y.genes),z.genes))))

outfile = paste("limma_edgeR_DESeq2_",args$contrast,"_vennDiagram.png",sep="")

png(outfile)
plot(vennD, doWeights = FALSE, type = "circles")
dev.off()

