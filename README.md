# R_samples
A small collection of R code I've written for various projects at CDC.

## Contents
Like the summary says, this is a very small selection of the work I've done in R for my statistical projects at CDC. Here's what I've uploaded so far:

1. `ga_pH.R` fits generalized estimating equations (GEEs) to some TB lab data, with the goal being to see whether the pH of gastric aspirate samples affects their positivity rate (i.e. yield) on MGIT culture. 
2. `mmwr.R` runs an interrupted time series (ITS) analysis on opioid use data. The data come from Blue Cross Blue Shield Massachusetts (BCBS MA), who implemented a policy change designed to affect opioid prescribing rates and wanted to see if it worked. The published paper is [here](https://www.cdc.gov/mmwr/volumes/65/wr/mm6541a1.htm).
3. `proportions.R` contains a bunch of custom functions for working with binomially-distributed data. It is poorly organized and leans heavily on Newcombe's score-based methods for constructing confidence intervals.
4. `sample_size.R` 

## Disclaimer
Most of the functions in these scripts should work, but some were abandoned during development and will thus crash, and others will work but fail to produce reasonable results.
