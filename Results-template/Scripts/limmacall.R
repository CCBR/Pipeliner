
## grab args
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
#FILE1 <- args[2]
#FILE2 <- args[3]
#CONTRASTS <- args[4]
# Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")
setwd(DIR) # new 
# rmarkdown::render("Scripts/LimmaReport.Rmd", params = list(
rmarkdown::render("LimmaReport.Rmd", params = list(
    folder = args[1],
    sampleinfo = args[2],
    data = args[3],
    contrasts = args[4],
    species = args[5],
    projectId = args[6],
    projectDesc = args[7],
    gtffile = args[8],
    dtype = args[9],
    karyobeds = args[10]
  ))

