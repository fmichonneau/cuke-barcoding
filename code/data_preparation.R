

source("code/functions_cukeBarcoding.R")
library(seqManagement)
library(doMC)
registerDoMC()
library(ape)
library(seqinr)
library(phylobase)


#### load CSV file
allDB <- read.csv(file="data/MARBoL_Echinos_VIII_2013.csv", stringsAsFactors=FALSE) # nrow = 7017
