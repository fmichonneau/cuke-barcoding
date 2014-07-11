
library(ape)
seqHol <- readRDS("data/cukeBarcodes-flagAmb.rds")
treH <- nj(dist.dna(seqHol, model="raw"))
saveRDS(treH, "data/cukeTree-raw.rds")
