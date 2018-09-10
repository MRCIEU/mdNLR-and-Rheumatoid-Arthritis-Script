###########################################
# mdNLR Rheumatoid arthritis
##########################################

source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(minfi)
library(GEOquery)
getGEOSuppFiles("GSE42861")
untar("GSE68777/GSE68777_RAW.tar", exdir = "GSE68777/idat")
head(list.files("GSE68777/idat", pattern = "idat"))

source("https://bioconductor.org/biocLite.R")
biocLite("minfi")
rgSet <- read.450k.exp("IDAT")


