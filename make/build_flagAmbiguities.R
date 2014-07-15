library(seqManagement)
library(ape)

algFile <- "data/cukeBarcodes-cleaned.fas"

### identify sequences with ambiguities and rename them
ambSeq <- checkAmbiguity(file=algFile)
oldNm <- names(ambSeq)
newNm <- paste(oldNm, "_", sapply(ambSeq, length), "amb", sep="")

seqHol <- read.dna(file=algFile, format="fasta")
dimnames(seqHol)[[1]][match(oldNm, dimnames(seqHol)[[1]])] <- newNm
seqHol <- cleanSeqLabels(seqHol, software="RAxML")

saveRDS(seqHol, file="data/cukeBarcodes-flagAmb.rds")
