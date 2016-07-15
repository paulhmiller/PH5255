# The 70-gene diagnostic signature for breast cancer progression risk
# genefu package

library(genefu)
data(sig.gene70)
dim(sig.gene70)
head(sig.gene70)[,1:6]

sum(is.na(sig.gene70$NCBI.gene.symbol))
grep("kinase", sig.gene70$Description)

# checking normality and transforming data
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

# exploring the genome
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
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
