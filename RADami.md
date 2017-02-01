```R
# RADami plot.locus.dist function

setwd()
library(RADami)

rubiset <- read.pyRAD("20170126_Rubi_27taxa.loci")
distancetable <- locus.dist(rubiset)

tree <- read.tree("RAxML_bipartitions.BestRapidBS_SNPs.tre")

plot.locus.dist(distancetable,tree, scalar = 5)

```