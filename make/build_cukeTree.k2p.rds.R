
args <- commandArgs(TRUE)
file.in <- args[1]
model <- args[2]
file.out <- args[3]

library(ape)
seqHol <- readRDS(file.in)
treH <- nj(dist.dna(seqHol, model=model))
bootH <- boot.phylo(treH, seqHol, function(xx) nj(dist.dna(xx), model=model), B=200)
saveRDS(treH, file.out)
