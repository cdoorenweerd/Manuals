# BAMM manual
By C. Doorenweerd. This manual assumes you have downloaded and compiled BAMM, that you have R (studio) and the packages APE, geiger, coda and BAMMtools installed (`install.packages("packagename")`). Load the packages using:

```R
library(ape)
library(geiger)
library(coda)
library(BAMMtools)
```

[BAMM CONFIGURATION FILE INFO](http://bamm-project.org/configuration.html)


# General useful R commands

Ensure you are working from the folder where your files are at:

```R
options(repos = c(CRAN = "http://cran.rstudio.com"))
getwd()
setwd("dir")
list.files()
rm(list = ls())
```

Simple tree handling commands:

```R
tree <- read.tree("treefile.tre")    # Reads a newick tree
outgroup = c("taxon1", "taxon2")    # Define outgroup taxa
prunedtree <- drop.tip(tree, outgroup)   # For removing outgroups
plot.phylo(tree, cex = 0.5)    # simple plot of tree
nodelabels(cex = 0.5, frame = "none")    # plots nodelabel numbers on tree
write.tree(MyTree, file="MyNewickTreefile.tre")    # Export a newick tree
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

## 0. Read in the BAMM output files

```R
tree <- read.tree("treefile.tre")    # Reads the newick tree
mcmcout <- read.csv("mcmc_out.txt", header = TRUE)   # Information on the MCMC run
edata <- getEventData(tree, eventdata = "bamm_event_data.txt", burnin=0.1)    # The main output from BAMM algorithm
```

## 1. Check if the run went well

Plot the log likelihood trace:
```R
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
plot.bammdata(edata, tau=0.001, breaksmethod='quantile', lwd=2, pal="temperature", legend=TRUE)
```

## 3. [Compute bayes factors](http://bamm-project.org/postprocess.html#bayes-factors-for-model-comparison)

To calculate the most likely number of rate shifts it is necessary to calculate bayes factors. To calculate M0 null model with simulation from prior only, in config file, set `simulatePriorShifts = 1`. Then compare the output from this run with the actual run, by calculating bayes factors:

```R

postfile <- "post_mcmc_out.txt"
bfmat <- computeBayesFactors(postfile, expectedNumberOfShifts=1, burnin=0.1)


postfile <- "BAMM_mcmc_out.txt"
priorfile <- "BAMM_prior_probs.txt"
bayesfactormatrix <- computeBayesFactors(postfile, priorfile, burnin = 0.1)
```

Differences of 20 are considered strong evidence for one scenario over another, differences of 50 as very strong.

## 4. Plotting individual scenario's

Select and plot the most likely scenario of number of shifts:

```R
summary(edata)
edata2 <- subsetEventData(edata, index = 2)
plot.bammdata(edata2, breaksmethod='quantile', lwd=2, pal="temperature", legend=TRUE)
addBAMMshifts(edata2, par.reset = FALSE, cex = 1)
```

## 5. [Branch specific shift probabilities](http://bamm-project.org/rateshifts.html#marginal-shift-probabilities)

To calculate marginal shift probabilities per branch:

```R
marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs, show.tip.label = FALSE)
```

However, these may be misleading. To convert these to bayes factor probabilities:

```R
branch_priors <- getBranchShiftPriors(tree, "BAMM_prior_probs.txt")
bf <- bayesFactorBranches(edata, branch_priors)
plot.phylo(bf, show.tip.label = FALSE)
```

Where the branch length is the bayes factor value. To see the values directly:

```R
bf$edge.length
```

## 6. [Bayesian 95% credible set of shift configurations](http://bamm-project.org/postprocess.html#bayesian-credible-sets-of-shift-configurations)


```R
priorshifts <- getBranchShiftPriors(tree, "m0prior_prior_probs.txt")    # The probs from the priors only run
css <- credibleShiftSet(edata, priorshifts, set.limit = 0.95)    # Set with 95% of credible shifts
summary(css)    # Lists the different scenarios
plot.credibleshiftset(css, lwd = 2, breaksmethod='quantile', pal = "temperature")    # Plots the different scenarios
```

## 7. [Evolutionary rate variation through time: grayscale](http://bamm-project.org/bammgraph.html#evolutionary-rate-variation-through-time-grayscale)

```R
starttime <- max(branching.times(tree))
plotRateThroughTime(edata, avgCol="black", start.time=starttime, ylim=c(0,1), cex.axis=2, intervalCol='gray80', intervals=c(0.05, 0.95), opacity=1)
```
Or specifying a single clade only:

```R
plotRateThroughTime(edata, avgCol="black", start.time=starttime, node=72, xlim=c(100,0), ylim=c(0,0.3), cex.axis=1, cex=1, intervalCol='gray80', intervals=c(0.05, 0.95), opacity=1)
```
