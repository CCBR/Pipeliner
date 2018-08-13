
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

