
source("R/parseRscriptArgs.R")
args <- commandArgs(TRUE)
args <- parseRscriptArgs(args)

file.in <- args$file.in
model <- args$model
file.out <- args$file.out
Nrep <- args$Nrep

set.seed(10101)

seqHol <- readRDS(file.in)
treH <- build_cukeTree(alg=seqHol, model=model, Nrep=Nrep)
saveRDS(treH, file.out)
