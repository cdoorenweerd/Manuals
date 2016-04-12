NOTE: This HOWTO assumes you have downloaded and compiled BAMM, that you have R (studio) and the packages APE, geiger, coda and BAMMtools installed (`install.packages("packagename")`). Load the packages using:

```R
library(ape)
library(geiger)
library(coda)
library(BAMMtools)
```

BAMM runs from R to prepare input files and analyse the results, the MCMC part runs from a C++ written program that has to be installed separately. The latter can be run on OpenStack and comes installed with the PhylOStack. All R parts are best run locally, using RStudio.


General useful R commands
=========================

```R
options(repos = c(CRAN = "http://cran.rstudio.com")) # sets repository
getwd() # shows working directory
setwd("dir") # sets working directory
list.files() # show files in working directory
rm(list = ls()) # clear loaded files
tree <- read.tree("treefile.tre")    # Reads a newick tree
outgroup = c("taxon1", "taxon2")    # Define outgroup taxa
prunedtree <- drop.tip(tree, outgroup)   # For removing outgroups
plot.phylo(tree, cex = 0.5)    # simple plot of tree
nodelabels(cex = 0.5, frame = "none")    # plots nodelabel numbers on tree
write.tree(MyTree, file="MyNewickTreefile.tre")    # Export a newick tree
```

Preparing input files
=====================

The MCMC part of BAMM requires two mandatory files and one optional:

- An ultrametric tree in NEWICK format, without outgroup taxa. Preferably with branch lengths that correspond with age in million of years. In the examples this file is named mytree.nwk and as an R object named as `tree`.

- a BAMM_configuration.txt file, see below

- optionally: a tab delimited text file defining the sampling percentage per clade, named samplingfractions.txt in the examples


An ultrametric newick tree
--------------------------

Although there are some options to create arbitrary ultrametric trees using the chronos function from APE in R, it is strongly recommended to start with a tree that has been created by BEAST.

1. Remove outgroups. Outgroup sampling is commonly far from complete, and in most cases it will be best to remove the outgroups before the analysis. This is easy to do in FigTree, where you select the ingroup clade, CTRL+C, open notepad and CTRL+V. Save the notepad file, it will contain the ingroup tree in nexus format. Open this file with FigTree and export the tree in newick format.

2. From R, load the tree and (double)check if the format is suitable for BAMM:

```R
library(ape)
tree <- read.tree("mytree.nwk")   # read your tree as tree object
is.ultrametric(tree)   # check for ultrametric
is.binary.tree(tree)   # check for binary
min(tree$edge.length)  # check all branchlengths are > 0
```

If all tests pass, your tree is ready. Define the tree name in BAMM_configuration.txt


Setting clade specific sampling percentages <http://bamm-project.org/advanced.html#accounting-for-non-random-incomplete-taxon-sampling-in-diversification-studies
-------------------------------------------

If your taxon names contain strange characters (dots, dashes etc), there will be single quotes around your taxon names, which need to be included in this taxonlist in order for BAMM and R to recognize the taxa. To get an accurate list of your taxa, use R:

```R
tree <- read.tree("mytree.nwk")
tiplist <- tree$tip.label
write.table(tiplist, "taxonlist.txt")
```

Open this list in Excel and separate the data using double quotes as separator, define the sampled percentages per clade as described in the manual and save as a tab-delimited text file. Make sure this text file has Linux style line breaks! 



BAMM_configuration.txt
----------------------

- First, use BAMMtools in R to determine the correct priors for your data:

```R
tree <- read.tree("mytree.nwk")
setBAMMpriors(tree)
```

This will create a text file in your working directory with values that you can copy paste into the BAMM_configuration.txt. You can use this template:

```
# BAMM configuration file for speciation/extinction analysis 
# ==========================================================
#
# Format
# ------
#
#     - Each option is specified as: option_name = option_value
#     - Comments start with # and go to the end of the line
#     - True is specified with "1" and False with "0"


################################################################################
# GENERAL SETUP AND DATA INPUT
################################################################################

modeltype = speciationextinction        
# Specify "speciationextinction" or "trait" analysis
                                  
treefile = finaltreekeep_noOG.nwk
# File name of the phylogenetic tree to be analyzed

runInfoFilename = run_info.txt
# File name to output general information about this run

sampleFromPriorOnly = 0                 
# Whether to perform analysis sampling from prior only (no likelihoods computed)

runMCMC = 1                             
# Whether to perform the MCMC simulation. If runMCMC = 0, the program will only
# check whether the data file can be read and the initial likelihood computed

simulatePriorShifts = 0
# Whether to simulate the prior distribution of the number of shift events,
# given the hyperprior on the Poisson rate parameter. This is necessary to
# compute Bayes factors

loadEventData = 0                       
# Whether to load a previous event data file

eventDataInfile = event_data_in.txt
# File name of the event data file to load, used only if loadEventData = 1

initializeModel = 1                     
# Whether to initialize (but not run) the MCMC. If initializeModel = 0, the
# program will only ensure that the data files (e.g., treefile) can be read

useGlobalSamplingProbability = 0        
# Whether to use a "global" sampling probability. If False (0), expects a file
# name for species-specific sampling probabilities (see sampleProbsFilename)
                                        
globalSamplingFraction = 0            
# The sampling probability. If useGlobalSamplingFraction = 0, this is ignored
# and BAMM looks for a file name with species-specific sampling fractions

sampleProbsFilename = samplingfractions.txt
# File name containing species-specific sampling fractions

# seed = 12345
# Seed for the random number generator. 
# If not specified (or is -1), a seed is obtained from the system clock

overwrite = 1
# If True (1), the program will overwrite any output files in the current
# directory (if present)


################################################################################
# PRIORS
################################################################################


###############################################

# Prior block chosen by BAMMtools::setBAMMpriors
# PASTE YOUR VALUES HERE

poissonRatePrior = 1.0

lambdaInitPrior = 3.90088932715874

lambdaShiftPrior = 0.0115198279154606

muInitPrior = 3.90088932715874

#### End Prior block
######################

lambdaIsTimeVariablePrior = 1
# Prior (probability) of the time mode being time-variable (vs. time-constant)


################################################################################
# MCMC SIMULATION SETTINGS & OUTPUT OPTIONS
################################################################################

numberOfGenerations = 2000000
# Number of generations to perform MCMC simulation

mcmcOutfile = mcmc_out.txt
# File name for the MCMC output, which only includes summary information about
# MCMC simulation (e.g., log-likelihoods, log-prior, number of processes)

mcmcWriteFreq = 1000
# Frequency in which to write the MCMC output to a file

eventDataOutfile = event_data.txt
# The raw event data (these are the main results). ALL of the results are
# contained in this file, and all branch-specific speciation rates, shift
# positions, marginal distributions etc can be reconstructed from this output.
# See R package BAMMtools for working with this output

eventDataWriteFreq = 1000
# Frequency in which to write the event data to a file

printFreq = 1000
# Frequency in which to print MCMC status to the screen

acceptanceResetFreq = 1000
# Frequency in which to reset the acceptance rate calculation
# The acceptance rate is output to both the MCMC data file and the screen

outName = BAMM
# Optional name that will be prefixed on all output files (separated with "_")
# If commented out, no prefix will be used


################################################################################
# OPERATORS: MCMC SCALING OPERATORS
################################################################################

updateLambdaInitScale = 2.0
# Scale parameter for updating the initial speciation rate for each process

updateLambdaShiftScale = 0.1
# Scale parameter for the exponential change parameter for speciation

updateMuInitScale = 2.0
# Scale parameter for updating initial extinction rate for each process

updateEventLocationScale = 0.05
# Scale parameter for updating LOCAL moves of events on the tree
# This defines the width of the sliding window proposal
 
updateEventRateScale = 4.0
# Scale parameter (proportional shrinking/expanding) for updating
# the rate parameter of the Poisson process 


################################################################################
# OPERATORS: MCMC MOVE FREQUENCIES
################################################################################

updateRateEventNumber = 0.1
# Relative frequency of MCMC moves that change the number of events

updateRateEventPosition = 1
# Relative frequency of MCMC moves that change the location of an event on the
# tree

updateRateEventRate = 1
# Relative frequency of MCMC moves that change the rate at which events occur 

updateRateLambda0 = 1
# Relative frequency of MCMC moves that change the initial speciation rate
# associated with an event

updateRateLambdaShift = 1
# Relative frequency of MCMC moves that change the exponential shift parameter
# of the speciation rate associated with an event

updateRateMu0 = 1
# Relative frequency of MCMC moves that change the extinction rate for a given
# event

updateRateLambdaTimeMode = 0
# Relative frequency of MCMC moves that flip the time mode
# (time-constant <=> time-variable)

localGlobalMoveRatio = 10.0
# Ratio of local to global moves of events 


################################################################################
# INITIAL PARAMETER VALUES
################################################################################

lambdaInit0 = 0.032
# Initial speciation rate (at the root of the tree)

lambdaShift0 = 0
# Initial shift parameter for the root process

muInit0 = 0.005
# Initial value of extinction (at the root)

initialNumberEvents = 0
# Initial number of non-root processes


################################################################################
# METROPOLIS COUPLED MCMC
################################################################################

numberOfChains = 4
# Number of Markov chains to run

deltaT = 0.01
# Temperature increment parameter. This value should be > 0
# The temperature for the i-th chain is calculated as 1 / [1 + deltaT * (i - 1)]

swapPeriod = 1000
# Number of generations in which to propose a chain swap

chainSwapFileName = chain_swap.txt
# File name in which to output data about each chain swap proposal.
# The format of each line is [generation],[rank_1],[rank_2],[swap_accepted]
# where [generation] is the generation in which the swap proposal was made,
# [rank_1] and [rank_2] are the chains that were chosen, and [swap_accepted] is
# whether the swap was made. The cold chain has a rank of 1.
 

################################################################################
# NUMERICAL AND OTHER PARAMETERS
################################################################################

minCladeSizeForShift = 1
# Allows you to constrain location of possible rate-change events to occur
# only on branches with at least this many descendant tips. A value of 1
# allows shifts to occur on all branches. 

segLength = 0.02
# Controls the "grain" of the likelihood calculations. Approximates the
# continuous-time change in diversification rates by breaking each branch into
# a constant-rate diversification segments, with each segment given a length
# determined by segLength. segLength is in units of the root-to-tip distance of
# the tree. So, if the segLength parameter is 0.01, and the crown age of your
# tree is 50, the "step size" of the constant rate approximation will be 0.5.
# If the value is greater than the branch length (e.g., you have a branch of
# length < 0.5 in the preceding example) BAMM will not break the branch into
# segments but use the mean rate across the entire branch.

```


Running BAMM
============

When you have prepared your input files you can run BAMM MCMC algorithm, which is as simple as typing the command below from the folder with the input files:

```sh
bamm -c BAMM_configuration.txt
```



Analysing the results with BAMMtools in R
=========================================

Load results
------------

1. Load the tree, if you haven't done so already

```R
tree <- read.tree("mytree.nwk")
```

2. Decide burnin percentage

```R
mcmcout <- read.csv("BAMM_mcmc_out.txt")   # Load the information on the MCMC run
plot(mcmcout$logLik ~ mcmcout$generation) # plot log likelihood vs generations
```

Typically your burnin will be somewhere between 0.1 and 0.25

3. Check if the post burnin effective sampling size is >>200

```R
burnstart <- floor(0.2 * nrow(mcmcout)) # choose a burnin here, currently 0.2
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
effectiveSize(postburn$N_shifts) # should be >>200
effectiveSize(postburn$logLik) # should be >>200
```

4. Load the main BAMM output: the event data

Use the determined burnin value.

```R
edata <- getEventData(tree, eventdata = "BAMM_event_data.txt", burnin=0.2)
summary(edata)
```


Visualize results
-----------------


- Mean phylorate plot (a quick look at the results)
```R
plot.bammdata(edata, tau=0.001, breaksmethod='quantile', lwd=2, pal="temperature", legend=TRUE)
```

- View 95% most credible shift set

```R
css <- credibleShiftSet(edata, expectedNumberOfShifts = 1, threshold=5, set.limit=0.95)
summary(css)
plot.credibleshiftset(css, plotmax=4)
```

- See how the shift probabilities are divided over clades

```R
marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs, show.tip.label = FALSE)
```
?this may be misleading


- Maximum a posteriori probability shift configuration

<bamm-project.org/rateshifts.html/#overall-best-shift-configuration>
Recommended for publication.

```R
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts = 1, threshold = 5)
plot.bammdata(best, lwd=1, breaksmethod='quantile', pal="temperature")
addBAMMshifts(best, cex=1)
```

- Evolutionary rate variation through time: 

    - Colour

```R
starttime <- max(branching.times(tree))
plotRateThroughTime(edata, intervalCol="darkgreen", avgCol="darkgreen", start.time=starttime, ylim=c(0,0.15))
```

    - grayscale 

<http://bamm-project.org/bammgraph.html#evolutionary-rate-variation-through-time-grayscale>

```R
starttime <- max(branching.times(tree))
plotRateThroughTime(edata, avgCol="black", start.time=starttime, ylim=c(0,1), cex.axis=2, intervalCol='gray80', intervals=c(0.05, 0.95), opacity=1)
```

    - Or specifying a single clade only:

```R
plotRateThroughTime(edata, avgCol="black", start.time=starttime, node=72, xlim=c(100,0), ylim=c(0,0.3), cex.axis=1, cex=1, intervalCol='gray80', intervals=c(0.05, 0.95), opacity=1)
```



######################
NOT SURE WHAT THESE DO
######################




## 3. [Compute bayes factors](http://bamm-project.org/postprocess.html#bayes-factors-for-model-comparison)

To calculate the most likely number of rate shifts it is necessary to calculate bayes factors. To calculate M0 null model with simulation from prior only, in config file, set `simulatePriorShifts = 1`. Then compare the output from this run with the actual run, by calculating bayes factors:

```R
postfile <- "BAMM_mcmc_out.txt"
priorfile <- "BAMM_prior_probs.txt"
bayesfactormatrix <- computeBayesFactors(postfile, priorfile, burnin = 0.1)
```

Differences of 20 are considered strong evidence for one scenario over another, differences of 50 as very strong.



- Plotting individual scenario's

Select and plot a specific scenario of number of shifts from the credible shift set

```R
summary(edata)
edata2 <- subsetEventData(edata, index = 2)
plot.bammdata(edata2, breaksmethod='quantile', lwd=2, pal="temperature", legend=TRUE)
addBAMMshifts(edata2, par.reset = FALSE, cex = 1)
```



## 6. [Bayesian 95% credible set of shift configurations](http://bamm-project.org/postprocess.html#bayesian-credible-sets-of-shift-configurations)


```R
priorshifts <- getBranchShiftPriors(tree, "m0prior_prior_probs.txt")    # The probs from the priors only run
css <- credibleShiftSet(edata, priorshifts, set.limit = 0.95)    # Set with 95% of credible shifts
summary(css)    # Lists the different scenarios
plot.credibleshiftset(css, lwd = 2, breaksmethod='quantile', pal = "temperature")    # Plots the different scenarios
```

rior_probs.txt")    # The probs from the priors only run
css <- credibleShiftSet(edata, priorshifts, set.limit = 0.95)    # Set with 95% of credible shifts
summary(css)    # Lists the different scenarios
plot.credibleshiftset(css, lwd = 2, breaksmethod='quantile', pal = "temperature")    # Plots the different scenarios
```

