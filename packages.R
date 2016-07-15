# Packages for Edx Bioconductor course

# As you work through exercises, if a Bioconductor package not found error is encountered, use
library(BiocInstaller)
biocLite([missing package name here, in quotes])


source("http://www.bioconductor.org/biocLite.R")
biocLite()
library(devtools)
install_github("genomicsclass/ph525x")

biocLite("bioBase")
#biocLite("GSE5859")
devtools::install_github("genomicsclass/GSE5859")
biocLite("annotate")
biocLite("BiocParallel")
biocLite("VariantAnnotation")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
biocLite("SNPlocs.Hsapiens.dbSNP.20120608")
biocLite("Homo.sapiens")



library(Biobase)
library(GSE5859)
library(annotate)
library(BiocParallel)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Homo.sapiens)







