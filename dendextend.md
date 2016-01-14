dendextend quick reference
==========================

Install denextend packages:

```R
options(repos = c(CRAN = "http://cran.rstudio.com"))
install.packages('dendextend')
install.packages('dendextendRcpp')
install.packages.2('colorspace')
```

Load required packages and set working dir:

```R
library(ape)
library("dendextend")
library("dendextendRcpp")
setwd("path/to/workingdir/")
```

Load newick trees. Trees must be rooted, no polytomies allowed.

```R
tree1 <- read.tree("firsttree.nwk")
tree2 <- read.tree("secondtree.nwk")
is.rooted(tree1)
is.rooted(tree2)
is.binary.tree(tree1)
is.binary.tree(tree2)
is.ultrametric(tree1)
is.ultrametric(tree2)
```

Make trees ultrametric (if not already so)

:code::R
```R
ctree1 <- chronos(tree1, lambda=1)
ctree2 <- chronos(tree2, lambda=1)
````

Compare the trees!

```R
tanglegram(ctree1,ctree2) # plots the two trees with colours
entanglement(ctree1,ctree2) # gives a numeric comparison value where lower is better
```
