Install and load libraries
==========================

```R
options(repos = c(CRAN = "http://cran.rstudio.com")) # sets repository
```

```R
install.packages("packagename") # install package
```

Load the packages using:

```R
library(ape)
library(geiger)
library(coda)
library(BAMMtools)
```


General useful R commands
=========================

```R
getwd() # shows working directory
setwd("dir") # sets working directory
list.files() # show files in working directory
rm(list = ls()) # clear loaded files
```

General tree handling commands
==============================

```R
tree <- read.tree("treefile.tre")    # Reads a newick tree
outgroup = c("taxon1", "taxon2")    # Define outgroup taxa
prunedtree <- drop.tip(tree, outgroup)   # For removing outgroups
plot.phylo(tree, cex = 0.5)    # simple plot of tree
nodelabels(cex = 0.5, frame = "none")    # plots nodelabel numbers on tree
write.tree(MyTree, file="MyNewickTreefile.tre")    # Export a newick tree
```

Create taxon list
=================

If your taxon names contain strange characters (dots, dashes etc), there will be single quotes around your taxon names, which need to be included in this taxonlist in order for BAMM and R to recognize the taxa. To get an accurate list of your taxa, use R:

```R
tree <- read.tree("mytree.nwk")
tiplist <- tree$tip.label
write.table(tiplist, "taxonlist.txt")
```

Open this list in Excel and separate the data using double quotes as separator, define the sampled percentages per clade as described in the manual and save as a tab-delimited text file. Make sure this text file has Linux style line breaks! 


Prune taxa from list
====================

Assuming you have a text file with a list of taxa that you want to remove from the tree, named :code:`droplist.txt`

```R
tree <- read.tree("mytree.nwk")

dropcsv <- read.csv("droplist.txt", header = FALSE)$V1
droplist <- as.character(dropcsv)
prunedtree <- drop.tip(tree, droplist)
```



