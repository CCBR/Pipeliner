library("rmarkdown")
args = commandArgs(trailingOnly=TRUE)
#arguments
#1=ChIPseeker.rmd ... fullpath
#2=bedFile
#3=genome ... hg19 or mm10
#4=output html file
#Sys.setenv(RSTUDIO_PANDOC="/usr/local/apps/Rstudio/0.98/bin/pandoc")
rmarkdown::render(args[1],
                  params = list(bedFile=args[2],genome=args[3]),
                  output_file=args[4])
