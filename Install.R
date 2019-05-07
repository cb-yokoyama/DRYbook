# インストール
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("GenomicRanges")

install.packages("devtools")
install.packages("changepoint")