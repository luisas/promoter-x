library(ips)
library(phangorn)

x <- read.fas("./data/promoters/sampled/proms_3.fa")
x <- as.phyDat(x)
y <- dist.hamming(x)
y
