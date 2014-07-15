
source("R/parseRscriptArgs.R")
args <- commandArgs(TRUE)
args <- parseRscriptArgs(args)

file.in <- args$file.in
model <- args$model
file.out <- args$file.out
Nrep <- args$Nrep

set.seed(10101)

library(ape)
seqHol <- readRDS(file.in)
treH <- nj(dist.dna(seqHol, model=model))
bootH <- boot.phylo(treH, seqHol, function(xx) nj(dist.dna(xx, model=model)), B=Nrep)
treH$node.label <- bootH
saveRDS(treH, file.out)
