library(changepoint)
source("R/cpttest.R")

# Import Methylation Data
H1_chr1 <- read.csv("extData/Lister2009/hg19/h19_H1_chr1.outfile")

# Domain demarcation in chr1
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
