
# MethylationExport
H1 <- read.csv("extData/Lister2009/hg19/h19_H1.outfile",sep="\t",head=F)
#IMR90 <- read.csv("h19_IMR90.outfile",sep="\t",head=F)


names(H1)<- c("seqid","pos","strand","mratio")
#names(IMR90)<- c("seqid","pos","strand","mratio")

Tmp <- split(H1, H1$seqid)
list2env(Tmp,env=.GlobalEnv)

H1_chr1 <- H1[which(H1$seqid=="chr1"),]

# > nrow(H1_chr20)
# [1] 656997
# > nrow(mc_h1_20f)
# [1] 933382
# > ncpts()
# 815
source("R/cpttest.R")
library(changepoint)


#chr1
cpt.H1_chr1.pv20 <- cpt.mean(H1_chr1$mratio,method="PELT",penalty = "Manual",pen.value = "2.0")
cpt.H1_chr1.pv18 <- cpt.mean(H1_chr1$mratio,method="PELT",penalty = "Manual",pen.value = "1.8")
cpt.H1_chr1.pv16 <- cpt.mean(H1_chr1$mratio,method="PELT",penalty = "Manual",pen.value = "1.6")
cpt.H1_chr1.pv14 <- cpt.mean(H1_chr1$mratio,method="PELT",penalty = "Manual",pen.value = "1.4")
cpt.H1_chr1.pv12 <- cpt.mean(H1_chr1$mratio,method="PELT",penalty = "Manual",pen.value = "1.2")
cpt.H1_chr1.pv10 <- cpt.mean(H1_chr1$mratio,method="PELT",penalty = "Manual",pen.value = "1.0")
cpt.H1_chr1.pv08 <- cpt.mean(H1_chr1$mratio,method="PELT",penalty = "Manual",pen.value = "0.8")
cpt.H1_chr1.pv06 <- cpt.mean(H1_chr1$mratio,method="PELT",penalty = "Manual",pen.value = "0.6")
cpt.H1_chr1.pv04 <- cpt.mean(H1_chr1$mratio,method="PELT",penalty = "Manual",pen.value = "0.4")
cpt.H1_chr1.pv02 <- cpt.mean(H1_chr1$mratio,method="PELT",penalty = "Manual",pen.value = "0.2")
