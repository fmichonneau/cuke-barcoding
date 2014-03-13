

source("code/functions_cukeBarcoding.R")
library(seqManagement)
library(doMC)
registerDoMC()
library(ape)
library(seqinr)
library(phylobase)

######## Converts XLSX spreadsheet into CSV
### See here for more info https://github.com/dagwieers/unoconv

system("unoconv -l&") # start listener
system("sleep 0.5;")
system("unoconv -f csv data/MARBoL_Echinos_VIII_2013.xlsx") # converts document
system("pkill unoconv") # kill process

#### load CSV file
allDB <- read.csv(file="data/MARBoL_Echinos_VIII_2013.csv", stringsAsFactors=FALSE) # nrow = 7017
