# Example Usgae: Rscript pcacall.R 'DEG_cntrl-test_0.5_2' 'STAR_files/sampletable.txt' 'DEG_cntrl-test_0.5_2/RawCountFile_RSEM_genes_filtered.txt' 'hg19seDEG' 'Enter CCBR Project Description and Notes here.'
## grab args
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
outHtml <- args[2]
# Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")
setwd(DIR) # new 
rmarkdown::render("PcaReport.Rmd",output_file=outHtml, params = list(
    folder = args[1],
    sampleinfo = args[3],
    data = args[4],
    projectId = args[5],
    projectDesc = args[6]
  ))

