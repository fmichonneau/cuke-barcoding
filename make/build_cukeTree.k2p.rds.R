
library(ape)
seqHol <- readRDS("data/cukeBarcodes-flagAmb.rds")
treH <- nj(dist.dna(seqHol, model="K80"))
saveRDS(treH, "data/cukeTree-k2p.rds")
