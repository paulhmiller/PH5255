# PH5255.5x Introduction to Bioconductor
# Week 1


# A review of what we measure
#============================

# The 70-gene diagnostic signature for breast cancer progression risk using genefu package
library(genefu)
data(sig.gene70)
dim(sig.gene70)
head(sig.gene70)[,1:6]
sum(is.na(sig.gene70$NCBI.gene.symbol))
grep("kinase", sig.gene70$Description)

# Checking for gaussian normality and transforming data
library(COPDSexualDimorphism.data)
data(lgrc.expr.meta)
sum(expr.meta$GENDER == "2-Female")
median(expr.meta$pkyrs)
summary(expr.meta$pkyrs)
shapiro.test(expr.meta$pkyrs)
qqnorm(expr.meta$pkyrs)
qqline(expr.meta$pkyrs)
boxplot(pkyrs~gender, data=expr.meta)
expr.meta$pyp1 = expr.meta$pkyrs+1
library(MASS)
lm1 = lm(pyp1~gender, data=expr.meta)
boxcox(lm1)
lambda <- 0.5
boxplot(I(pyp1^lambda)~gender, data=expr.meta)

# Exploring the genome
library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome.Hsapiens.UCSC.hg19
chr11seq <- BSgenome.Hsapiens.UCSC.hg19[["chr11"]]
subseq(chr11seq,start=10^6,width=25)
?countPattern
countPattern("ATG", chr11seq)
countPattern("TGA", chr11seq)
countPattern("TAA", chr11seq)
countPattern("TAG", chr11seq)
?alphabetFrequency
chr7seq <- BSgenome.Hsapiens.UCSC.hg19[["chr7"]]
alphabetFrequency(chr7seq)
31671670 / sum(alphabetFrequency(chr7seq))

# Locations of SNPs in humans 
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
snps144 <- SNPlocs.Hsapiens.dbSNP144.GRCh37
s17 <- snplocs(snps144, "ch17")
head(s17)
#What is the location on chr17 of SNP rs73971683?
s17[s17$RefSNP_id == "73971683", ]

# GWAS: Linking SNP genotypes to disease risk 
library(gwascat)
data(ebicat37)
ebicat37
sort(table(ebicat37$CHR_ID),decreasing=TRUE)
head(sort(table(ebicat37$DISEASE.TRAIT),decreasing=TRUE))

# Gene expression with microarrays
library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e[,1:5])
table(tissue)
boxplot(e["209169_at", ] ~ tissue, las=2)
IDs = c("201884_at", "209169_at", "206269_at", "207437_at", "219832_s_at", "212827_at")
par(mfrow=c(2,3))
for (i in IDs){
  boxplot(e[i, ] ~ tissue, las=2, main=i)
}

# Discovery of microarray annotation in Bioconductor, search for chip ID on bioconductor
# Oligo sequences on affymetrix arrays 
library(BiocInstaller)
biocLite("hgu133aprobe")
library(hgu133aprobe)
head(hgu133aprobe)
nrow(hgu133aprobe[hgu133aprobe$Probe.Set.Name=="206269_at",])

