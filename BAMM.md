# BAMM manual
By C. Doorenweerd. This manual assumes you have downloaded and compiled BAMM, that you have R (studio) and the packages APE, coda and BAMMtools installed.

[BAMM CONFIGURATION FILE INFO](http://bamm-project.org/configuration.html)


# General useful R commands

```R
options(repos = c(CRAN = "http://cran.rstudio.com"))
getwd()
setwd("dir")
list.files()
rm(list = ls())

drop.tip(tree, "tipname")   # For removing outgroups
cladeA <- node.leaves(tree, mrca(tree)["taxon1", "taxon2"])
plot.phylo(tree)
write.tree(MyTree, file="MyNewickTreefile.tre")    # Will export a newick tree
```

# Running BAMM

## A: Speciation-extinction analyses

Check your tree is in the correct format using the APE package in R:
```R
library(ape)
tree <- read.tree("mytree.tre")   # read your tree as tree object
is.ultrametric(tree)   # check for ultrametric
is.binary.tree(tree)   # check for binary
min(tree$edge.length)  # check all branchlengths are > 0
```

Choose the priors using BAMMtools:
```R
library(BAMMtools)
setBAMMpriors(tree)
```
and paste the values from the created priors.txt into the configuration.txt


To account for random incomplete taxon sampling in diversification studies, make sure to set in the config file:

``` 
useGlobalSamplingProbability = 1
globalSamplingFraction = 0.50    # assumes 50% of all species sampled
```

Set the input tree file name in the config file
set the number of gens etc. `eventDataWriteFreq` value should equate to 1000 - 5000 output samples. Make sure `simulatePriorShifts = 0` for the actual run.

Run BAMM with:
```sh
bamm -c configuration.txt
```

## [B: Phenotypic evolution](http://bamm-project.org/quickstart.html#phenotypic-evolution) 

Much like A, except in the configfile `modeltype = trait`. And it needs a tab delimited trait file with taxonnames in first column and character states in the second. Define the name of this trait file in the config file.


## [C: MEDUSA-like model](http://bamm-project.org/advanced.html#medusa-like-model)

Rather than comparing clades, it looks for general increases or slow-downs in diversification rates (I THINK). To start this, use these settings in the configuration file:

```
updateRateEventNumber = 0.1
updateRateEventPosition = 1
updateRateEventRate = 1
updateRateLambda0 = 1
updateRateLambdaShift = 0
updateRateMu0 = 1
lambdaShift0 = 0
```


# BAMMtools workflow

## 1. Check if the run went well

In R, Plot the log likelihood trace:
```R
mcmcout <- read.csv("mcmc_out.txt", header = True)
plot(mcmcout$logLik ~ mcmcout$generation)
```

Discard burnin and examine ESS values and shift acceptance:
```R
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
```

## 2. Mean phylorate plot (a quick look at the results)
```R
library(BAMMtools)
tree <- read.tree("mytree.tre")
edata <- getEventData(tree, eventdata = "bamm_event_data.txt", burnin=0.1)
plot.bammdata(edata, tau=0.001, breaksmethod='quantile', lwd=2, pal="temperature", legend=TRUE)
```

## 3. [How many distinct rate shifts?](http://bamm-project.org/postprocess.html#how-many-rate-shifts)

Create a bammdata object called "edata" from our eventdata and summarize rateshifts:
```R
library(BAMMtools)
tree <- read.tree("mytree.tre")
edata <- getEventData(tree, eventdata = "bamm_event_data.txt", burnin=0.1)
summary(edata)
```

## 4. [Compute bayes factors](http://bamm-project.org/postprocess.html#bayes-factors-for-model-comparison)

Next to the run, calculate M0 null model with simulation from prior only. To get this, in config file set `simulatePriorShifts = 1` and `outName = M0prior`. Then compare the output from this run with the actual run, by calculating bayes factors:

```R
postfile <- "BAMM_mcmc_out.txt"
priorfile <- "M0prior_mcmc_out.txt"
bayesfactormatrix <- computeBayesFactors(postfile, priorfile, burnin = 0.1)
```

## 5. [Bayesian 95% credible set of shift configurations](http://bamm-project.org/postprocess.html#bayesian-credible-sets-of-shift-configurations)


```R
library(BAMMtools)
tree <- read.tree("mytree.tre")
edata <- getEventData(tree, eventdata = "bamm_event_data.txt", burnin=0.1)

priorshifts <- getBranchShiftPriors(tree, "m0prior_prior_probs.txt")
css <- credibleShiftSet(edata, priorshifts, set.limit = 0.95)
summary(css)
plot.credibleshiftset(css, lwd = 2, breaksmethod='quantile', pal = "temperature")
```

## 6. [Evolutionary rate variation through time: grayscale](http://bamm-project.org/bammgraph.html#evolutionary-rate-variation-through-time-grayscale)

```R
plot.new()
tree <- read.tree("mytree.tre")
starttime <- max(branching.times(tree))
edata <- getEventdata(tree, "BAMM_event_data.txt", burnin = 0.1)
plotRateThroughTime(edata, avgCol="black", start.time=starttime, ylim=c(0,1), cex.axis=2, intervalCol='gray80', intervals=c(0.05, 0.95), opacity=1)
```



